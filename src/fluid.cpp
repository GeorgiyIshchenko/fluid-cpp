#include "field_reader.hpp"
#include "types.hpp"
#include <array>
#include <barrier>
#include <bits/stdc++.h>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

template <typename vFlowType> struct BFSNode
{
    int x, y;
    vFlowType lim;
};

size_t global_n;
size_t global_m;

template <typename pType, typename vType, typename vFlowType, int N, int M>
class FluidSim
{
public:
    size_t n;
    size_t m;
    bool simProp = false;
    int num_threads = 1;
    bool alive = true;

    std::vector<thread> workers;
    std::barrier<std::function<void(void)>> start_barrier, end_barrier;

    PropagateFlowBorders borders;

    constexpr static bool sizeMatch = (N * M) != 0;

    constexpr static size_t T = 1'000'000;
    constexpr inline static std::array<pair<int, int>, 4> deltas = {
        { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } }
    };

    using fieldType = conditional_t<sizeMatch, array<array<char, M>, N>,
                                    vector<vector<char>>>;

    inline static fieldType field{};

    inline static vType rho[256]{};

    inline static conditional_t<sizeMatch, array<array<pType, M>, N>,
                                vector<vector<pType>>>
        p{}, old_p{};

    template <typename T> struct VectorField
    {
        conditional_t<sizeMatch, array<array<array<T, deltas.size()>, M>, N>,
                      vector<vector<vector<T>>>>
            v{};

        VectorField()
            requires(!sizeMatch)
        {
            v.resize(global_n);
            for (size_t i = 0; i < global_n; ++i)
            {
                v[i].resize(global_m);
                for (size_t j = 0; j < global_m; ++j)
                {
                    v[i][j].resize(deltas.size());
                }
            }
        }

        VectorField()
            requires(sizeMatch)
        = default;

        template <typename V> T& add(int x, int y, int dx, int dy, V dv)
        {
            return get(x, y, dx, dy) += dv;
        }

        T& get(int x, int y, int dx, int dy)
        {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    inline static VectorField<vType> velocity{}, velocity_flow{};
    inline static conditional_t<sizeMatch, array<array<int, M>, N>,
                                vector<vector<int>>>
        last_use{};
    inline static int UT = 0;

    conditional_t<sizeMatch, array<array<pType, M>, N>, vector<vector<pType>>>
        dirs{};

    mt19937 rnd{ 1337 };

    void init_workers()
    {
        for (auto&& border : borders.borders)
        {
            workers.emplace_back(
                [this, &border]()
                {
                    while (alive)
                    {
                        start_barrier.arrive_and_wait();
                        // cout << "WORKER" << "\n";
                        for (size_t x = 0; x < n; ++x)
                        {
                            for (size_t y = (size_t)border.first; y < (size_t)border.second; ++y)
                            {
                                if (field[x][y] != '#' && last_use[x][y] != UT)
                                {
                                    auto [t, local_prop, _] =
                                        propagate_flow(x, y, 1, true, border);
                                    if (t > 0.)
                                    {
                                        simProp = 1;
                                    }
                                }
                            }
                        }
                        end_barrier.arrive_and_wait();
                    }
                });
        }
    }

    FluidSim(FieldReader& inpField, int num_threads)
        requires(sizeMatch)
        : num_threads(num_threads),
          start_barrier(num_threads + 1, []() {}),
          end_barrier(num_threads + 1, []() {}),
          borders(inpField.M, num_threads)
    {
        n = inpField.N;
        m = inpField.M;
        cout << "STATIC FLUID CTOR: " << n << " " << m << " "
             << inpField.field.size() << " " << inpField.field[0].size()
             << endl;
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                cout << inpField.field[i][j];
            }
            cout << endl;
        }
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                field[i][j] = inpField.field[i][j];
            }
        }

        init_workers();

        cout << "CTOR SUCCESSED!" << endl;
    }

    FluidSim(FieldReader& inpField, int num_threads)
        requires(!sizeMatch)
        : num_threads(num_threads),
          start_barrier(num_threads, []() {}),
          end_barrier(num_threads, []() {}),
          borders(inpField.M, num_threads)
    {
        n = inpField.N;
        m = inpField.M;
        cout << "DYNAMIC FLUID CTOR: " << n << " " << m << " "
             << inpField.field.size() << " " << inpField.field[0].size()
             << endl;
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                cout << inpField.field[i][j];
            }
            cout << endl;
        }
        field.resize(n, vector<char>(m));
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                field[i][j] = inpField.field[i][j];
            }
        }

        init_workers();

        p.resize(n, vector<pType>(m, 0));
        old_p.resize(n, vector<pType>(m, 0));
        last_use.resize(n, vector<int>(m, 0));
        dirs.resize(n, vector<pType>(m, 0));
    }

    ~FluidSim()
    {
        alive = false;
    }

    tuple<vFlowType, bool, pair<int, int>>
    propagate_flow(int x, int y, vFlowType lim, bool checkBorders = false,
                   PropagateFlowBorder border = { 0, 0 })
    {
        last_use[x][y] = UT - 1;
        vFlowType ret = 0;
        for (auto [dx, dy] : deltas)
        {
            int nx = x + dx, ny = y + dy;

            if (checkBorders)
            {
                if (ny < border.first ||
                    ny >= border.second)
                {
                    continue;
                }
            }

            if (field[nx][ny] != '#' && last_use[nx][ny] < UT)
            {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap)
                {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(lim, static_cast<vFlowType>(cap - flow));
                if (last_use[nx][ny] == UT - 1)
                {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    // cerr << "A " << x << " " << y << " -> " << nx << " " <<
                    // ny
                    //      << " " << vp << " / " << lim << "\n";
                    return { vp, 1, { nx, ny } };
                }
                auto [t, prop, end] =
                    propagate_flow(nx, ny, vp, checkBorders, border);
                ret += t;
                if (prop)
                {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use[x][y] = UT;
                    // cerr << "B " << x << " " << y << " -> " << nx << " " <<
                    // ny
                    //      << " " << t << " / " << lim << "\n";
                    return { t, prop && end != pair(x, y), end };
                }
            }
        }
        last_use[x][y] = UT;
        return { ret, 0, { 0, 0 } };
    }

    double random01()
    {
        static std::uniform_real_distribution<double> distr(0, 1);
        return distr(rnd);
    }

    void propagate_stop(int x, int y, bool force = false)
    {
        if (!force)
        {
            bool stop = true;
            for (auto [dx, dy] : deltas)
            {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
                    velocity.get(x, y, dx, dy) > 0.)
                {
                    stop = false;
                    break;
                }
            }
            if (!stop)
            {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas)
        {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT ||
                velocity.get(x, y, dx, dy) > 0.)
            {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    vType move_prob(int x, int y)
    {
        vType sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i)
        {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT)
            {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0.)
            {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    struct ParticleParams
    {
        char type;
        pType cur_p;
        // array<vType, deltas.size()> v;

        inline static conditional_t<sizeMatch, array<vType, deltas.size()>,
                                    vector<vType>>
            v;

        ParticleParams()
            requires(sizeMatch)
        = default;

        ParticleParams()
            requires(!sizeMatch)
        {
            v.resize(deltas.size());
        }

        void swap_with(int x, int y)
        {
            swap(field[x][y], type);
            swap(p[x][y], cur_p);
            swap(velocity.v[x][y], v);
        }
    };

    bool propagate_move(int x, int y, bool is_first)
    {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do
        {
            std::array<pType, deltas.size()> tres;
            pType sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i)
            {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT)
                {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0.)
                {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0.)
            {
                break;
            }

            pType p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0. && field[nx][ny] != '#' &&
                   last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i)
        {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
                velocity.get(x, y, dx, dy) < 0.)
            {
                propagate_stop(nx, ny);
            }
        }
        if (ret)
        {
            if (!is_first)
            {
                ParticleParams pp{};
                pp.swap_with(x, y);
                pp.swap_with(nx, ny);
                pp.swap_with(x, y);
            }
        }
        return ret;
    }
};

template <typename pType, typename vType, typename vFlowType, size_t N,
          size_t M>
void do_main(FluidSim<pType, vType, vFlowType, N, M>& sim)
{
    cout << "Starting simultaion..." << endl;
#ifdef BENCH
    auto start = chrono::high_resolution_clock::now();
#endif
    sim.rho[' '] = 0.01;
    sim.rho['.'] = 1000;
    pType g = 0.1;

    for (size_t x = 0; x < sim.n; ++x)
    {
        for (size_t y = 0; y < sim.m; ++y)
        {
            if (sim.field[x][y] == '#')
                continue;
            for (auto [dx, dy] : sim.deltas)
            {
                sim.dirs[x][y] +=
                    static_cast<pType>(sim.field[x + dx][y + dy] != '#');
            }
        }
    }

    for (size_t i = 0; i < sim.T; ++i)
    {

        pType total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < sim.n; ++x)
        {
            for (size_t y = 0; y < sim.m; ++y)
            {
                if (sim.field[x][y] == '#')
                    continue;
                if (sim.field[x + 1][y] != '#')
                    sim.velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        // memcpy(sim.old_p, sim.p, sizeof(sim.p));
        for (size_t i = 0; i < sim.n; ++i)
        {
            sim.old_p[i] = sim.p[i];
        }

        for (size_t x = 0; x < sim.n; ++x)
        {
            for (size_t y = 0; y < sim.m; ++y)
            {
                if (sim.field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : sim.deltas)
                {
                    int nx = x + dx, ny = y + dy;
                    if (sim.field[nx][ny] != '#' &&
                        sim.old_p[nx][ny] < sim.old_p[x][y])
                    {
                        auto delta_p = sim.old_p[x][y] - sim.old_p[nx][ny];
                        auto force = delta_p;
                        auto& contr = sim.velocity.get(nx, ny, -dx, -dy);
                        if (contr * sim.rho[(int)sim.field[nx][ny]] >= force)
                        {
                            contr -= force / sim.rho[(int)sim.field[nx][ny]];
                            continue;
                        }
                        force -= static_cast<pType>(
                            contr * sim.rho[(int)sim.field[nx][ny]]);
                        contr = 0;
                        sim.velocity.add(x, y, dx, dy,
                                         force / sim.rho[(int)sim.field[x][y]]);
                        sim.p[x][y] -= force / sim.dirs[x][y];
                        total_delta_p -= force / sim.dirs[x][y];
                    }
                }
            }
        }

        // Make flow from velocities
        sim.velocity_flow = {};
        sim.simProp = false;
        do
        {
            sim.UT += 2;

            sim.simProp = false;

            sim.start_barrier.arrive_and_wait();
            sim.end_barrier.arrive_and_wait();

            for (auto&& border : sim.borders.borders)
            {
                for (size_t x = 0; x < sim.n; ++x)
                {
                    if ((size_t)border.second < sim.m)
                    {
                        if (sim.field[x][border.second] != '#' &&
                            sim.last_use[x][border.second] != sim.UT)
                        {
                            auto [t, local_prop, _] =
                                sim.propagate_flow(x, border.second, 1, false);
                            if (t > 0.)
                            {
                                sim.simProp = 1;
                            }
                        }
                    }
                }
            }

        } while (sim.simProp);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < sim.n; ++x)
        {
            for (size_t y = 0; y < sim.m; ++y)
            {
                if (sim.field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : sim.deltas)
                {
                    auto old_v = sim.velocity.get(x, y, dx, dy);
                    auto new_v = sim.velocity_flow.get(x, y, dx, dy);
                    if (old_v > 0.)
                    {
                        assert(new_v <= old_v);
                        sim.velocity.get(x, y, dx, dy) = new_v;
                        auto force =
                            (old_v - new_v) * sim.rho[(int)sim.field[x][y]];
                        if (sim.field[x][y] == '.')
                            force *= 0.8;
                        if (sim.field[x + dx][y + dy] == '#')
                        {
                            sim.p[x][y] += force / sim.dirs[x][y];
                            total_delta_p += force / sim.dirs[x][y];
                        }
                        else
                        {
                            sim.p[x + dx][y + dy] +=
                                force / sim.dirs[x + dx][y + dy];
                            total_delta_p += force / sim.dirs[x + dx][y + dy];
                        }
                    }
                }
            }
        }

        sim.UT += 2;
        sim.simProp = false;
        for (size_t x = 0; x < sim.n; ++x)
        {
            for (size_t y = 0; y < sim.m; ++y)
            {
                if (sim.field[x][y] != '#' && sim.last_use[x][y] != sim.UT)
                {
                    if (sim.random01() < sim.move_prob(x, y))
                    {
                        sim.simProp = true;
                        sim.propagate_move(x, y, true);
                    }
                    else
                    {
                        sim.propagate_stop(x, y, true);
                    }
                }
            }
        }

        if (sim.simProp)
        {
            cout << "Tick " << i << ":\n";
            for (size_t x = 0; x < sim.n; ++x)
            {
                for (size_t y = 0; y < sim.m; ++y)
                {
                    cout << sim.field[x][y];
                }
                cout << "\n";
            }
        }

#ifdef BENCH
        if (i >= 1000)
        {
            auto now = chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = now - start;
            cout << "Time elapsed: " << elapsed.count() << " s\n";            break;
        }
#endif
    }
}

void printDelim()
{
    cout << "#############################################" << endl;
}

int main(int argc, char** argv)
{

    // Outpud types debug info

    printDelim();

    cout << "Types string: " << TYPES_STRING << endl;

    cout << "Types amount: " << TYPES_COUNT << endl;

    cout << "Type list: ";

    for (auto&& s : typesStrArr)
    {
        cout << s << " ";
    }

    cout << endl;

    // Outpud sizes debug info

    printDelim();

    cout << "Sizes string: " << SIZES_STRING << endl;

    cout << "Sizes amount: " << SIZES_COUNT << endl;

    // Parse args

    string pTypeStr;
    string vTypeStr;
    string vFlowTypeStr;
    string fieldPath;
    int numThreads;
    if (argc < 4)
    {
        cout << "Not enogh arguments" << endl;
        exit(1);
    }
    else
    {
        for (int i = 0; i < argc; ++i)
        {
            string arg{ argv[i] };
            if (arg.substr(0, 2) == "--")
            {
                string paramType = {};
                size_t arg_idx = 2;
                while (arg[arg_idx] != '=')
                {
                    paramType += arg[arg_idx];
                    ++arg_idx;
                }
                arg_idx++;

                string value;
                while (arg_idx < arg.size())
                {
                    value += arg[arg_idx];
                    ++arg_idx;
                }

                if (paramType == "field")
                {
                    fieldPath = value;
                    continue;
                }

                if (paramType == "j")
                {
                    numThreads = stoi(value);
                    continue;
                }

                value = strToTypeStr(value);

                if (paramType == "p-type")
                {
                    pTypeStr = value;
                }
                else if (paramType == "v-type")
                {
                    vTypeStr = value;
                }
                else if (paramType == "v-flow-type")
                {
                    vFlowTypeStr = value;
                }
            }
        }
    }

    printDelim();

    cout << "Types that will be used: " << endl;
    cout << "p-type: " << pTypeStr << endl;
    cout << "v-type: " << vTypeStr << endl;
    cout << "v-flow-type: " << vFlowTypeStr << endl;

    assert(!pTypeStr.empty());
    assert(!vTypeStr.empty());
    assert(!vFlowTypeStr.empty());

    // Create field

    printDelim();

    FieldReader fr{ fieldPath };

    global_n = fr.N;
    global_m = fr.M;

    cout << "Field read successfully" << endl;
    cout << "Field sizes: " << fr.N << "x" << fr.M << endl;

    // Thread num

    printDelim();

    cout << "Number of threads: " << numThreads << endl;

    // Choose instance

    printDelim();

    bool sizeMatches = false;

    strToType(
        pTypeStr,
        [&]<typename pType>
        {
            cout << "pType chosen successfully: " << pTypeStr << endl;
            strToType(
                vTypeStr,
                [&]<typename vType>
                {
                    cout << "vType chosen successfully: " << vTypeStr << endl;
                    strToType(
                        vFlowTypeStr,
                        [&]<typename vFlowType>
                        {
                            cout << "vFlowType chosen successfully: "
                                 << vFlowTypeStr << endl;
                            printDelim();
                            numsToSize(fr.N, fr.M,
                                       [&]<typename sizeType>
                                       {
                                           sizeMatches = true;
                                           cout << "Sizes matched: " << fr.N
                                                << "x" << fr.M << endl;
                                           FluidSim<pType, vType, vFlowType,
                                                    sizeType::n, sizeType::m>
                                               sim(fr, numThreads);
                                           do_main<pType, vType, vFlowType,
                                                   sizeType::n, sizeType::m>(
                                               sim);
                                       });

                            if (!sizeMatches)
                            {
                                cout << "Sizes do not match" << endl;
                                FluidSim<pType, vType, vFlowType, 0, 0> sim(
                                    fr, numThreads);
                                do_main<pType, vType, vFlowType, 0, 0>(sim);
                            }
                        });
                });
        });

    return 0;
}
