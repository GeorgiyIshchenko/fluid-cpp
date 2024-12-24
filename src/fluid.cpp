#include "adapting_array.hpp"
#include "field_reader.hpp"
#include "types.hpp"
#include <array>
#include <atomic>
#include <barrier>
#include <bits/stdc++.h>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iostream>
#include <omp.h>
#include <string>
#include <string_view>
#include <thread>
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
    // General

    size_t n;
    size_t m;

    constexpr static bool sizeMatch = (N * M) != 0;

    size_t T = 1'000'000;
    int UT = 0;
    constexpr static std::array<pair<int, int>, 4> deltas = {
        { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } }
    };

    // Threads

    int num_threads = 1;
    atomic<bool> prop = false;
    bool alive = true;
    std::vector<jthread> workers;
    std::barrier<std::function<void(void)>> start_barrier, end_barrier;

    PropagateFlowBorders borders;

    // Arrays

    using fieldType = AdaptingArray<char, N, M, sizeMatch>;

    fieldType field{};

    vType rho[256]{};

    AdaptingArray<pType, N, M, sizeMatch> p, old_p;

    using last_uded_t = AdaptingArray<int, N, M, sizeMatch>;
    last_uded_t last_use, last_use_copy;

    AdaptingArray<pType, N, M, sizeMatch> dirs;

    template <typename T> struct VectorField
    {
        AdaptingArray<array<T, deltas.size()>, N, M, sizeMatch> v;

        VectorField()
            requires(!sizeMatch)
            : v(global_n, global_m, { 0 })
        {
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
            return v(x, y)[i];
        }
    };

    VectorField<vType> velocity, velocity_flow;

    mt19937 rnd{ 1337 };

    // DEBUG

    vector<int> debug_thread;

    vector<vector<char>> debug_thread_matrix = {};

    std::chrono::duration<double> debug_elapsed_waiting = {};

    void init_workers()
    {
        debug_thread.resize(borders.borders.size());
        debug_thread_matrix.resize(n, vector<char>(m, '-'));
        for (size_t x = 0; x < n; ++x)
        {
            for (size_t y = 0; y < m; ++y)
            {
                cout << debug_thread_matrix[x][y];
            }
            cout << "\n";
        }
        int i = 0;
        for (auto border : borders.borders)
        {
            workers.emplace_back(
                [this, border, i]()
                {
                    int c = i;
                    while (alive)
                    {
                        auto st = chrono::high_resolution_clock::now();
                        start_barrier.arrive_and_wait();
                        debug_elapsed_waiting +=
                            chrono::high_resolution_clock::now() - st;
                        for (int x = 0; x < static_cast<int>(n); ++x)
                        {
                            for (int y = border.first; y < border.second; ++y)
                            {
                                if (field(x, y) != '#' &&
                                    last_use_copy(x, y) != UT)
                                {
                                    auto [t, local_prop, _] = propagate_flow(
                                        x, y, 1, last_use_copy, true, border);
                                    if (t > 0.)
                                    {
                                        debug_thread[c]++;
                                        debug_thread_matrix[x][y] =
                                            to_string(c)[0];
                                        prop = true;
                                    }
                                }
                            }
                        }
                        st = chrono::high_resolution_clock::now();
                        end_barrier.arrive_and_wait();
                        debug_elapsed_waiting +=
                            chrono::high_resolution_clock::now() - st;
                    }
                });
            i++;
        }
    }

    FluidSim(FieldReader& inpField, int num_threads)
        requires(sizeMatch)
        : n(inpField.N),
          m(inpField.M),
          num_threads(num_threads),
          start_barrier(num_threads + 1, []() {}),
          end_barrier(num_threads + 1, []() {}),
          borders(inpField.M, num_threads),
          field(),
          p(),
          old_p(),
          last_use(),
          last_use_copy(),
          dirs()
    {
        n = inpField.N;
        m = inpField.M;

        cout << "STATIC FLUID CTOR: " << n << " " << m << " "
             << inpField.field.size() << " " << inpField.field[0].size()
             << "\n";
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                cout << inpField.field[i][j];
            }
            cout << "\n";
        }
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                field(i, j) = inpField.field[i][j];
            }
        }

        init_workers();

        cout << "CTOR SUCCESSED!" << "\n";
    }

    FluidSim(FieldReader& inpField, int num_threads)
        requires(!sizeMatch)
        : n(inpField.N),
          m(inpField.M),
          num_threads(num_threads),
          start_barrier(num_threads + 1, []() {}),
          end_barrier(num_threads + 1, []() {}),
          borders(inpField.M, num_threads),
          field(n, m),
          p(n, m),
          old_p(n, m),
          last_use(n, m),
          last_use_copy(n, m),
          dirs(n, m)
    {
        cout << "DYNAMIC FLUID CTOR: " << n << " " << m << " "
             << inpField.field.size() << " " << inpField.field[0].size()
             << "\n";
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                cout << inpField.field[i][j];
            }
            cout << "\n";
        }
        field = { n, m };
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                field(i, j) = inpField.field[i][j];
            }
        }

        init_workers();

        cout << "CTOR SUCCESSED!" << "\n";
    }

    FluidSim(const FluidSim& other) = delete;
    FluidSim(FluidSim&& other) = delete;

    pair<size_t, size_t> debug_missed = { 0, 0 };

    ~FluidSim()
    {
        alive = false;
    }

    tuple<vFlowType, bool, pair<int, int>>
    propagate_flow(int x, int y, vFlowType lim, last_uded_t& last_use,
                   bool checkBorders = false,
                   PropagateFlowBorder border = { 0, 0 })
    {

        if (checkBorders)
        {
            if (y < border.first || y >= border.second)
            {
                debug_missed.first++;
                return {};
            }
        }

        last_use(x, y) = UT - 1;
        vFlowType ret = 0;
        for (auto [dx, dy] : deltas)
        {
            int nx = x + dx, ny = y + dy;

            if (checkBorders)
            {
                if (ny < border.first || ny >= border.second)
                {
                    debug_missed.first++;
                    continue;
                }
            }
            debug_missed.second++;

            if (field(nx, ny) != '#' && last_use(nx, ny) < UT)
            {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap)
                {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(lim, static_cast<vFlowType>(cap - flow));
                if (last_use(nx, ny) == UT - 1)
                {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use(x, y) = UT;
                    // cerr << "A " << x << " " << y << " -> " << nx << " " <<
                    // ny
                    //      << " " << vp << " / " << lim << "\n";
                    return { vp, 1, { nx, ny } };
                }
                auto [t, prop, end] =
                    propagate_flow(nx, ny, vp, last_use, checkBorders, border);
                ret += t;
                if (prop)
                {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use(x, y) = UT;
                    // cerr << "B " << x << " " << y << " -> " << nx << " " <<
                    // ny
                    //      << " " << t << " / " << lim << "\n";
                    return { t, prop && end != pair(x, y), end };
                }
            }
        }
        last_use(x, y) = UT;
        return { ret, 0, { 0, 0 } };
    }

    double random01()
    {
        static std::uniform_real_distribution<double> distr(0, 1);
        return 1 - distr(rnd);
    }

    void propagate_stop(int x, int y, bool force = false)
    {
        if (!force)
        {
            bool stop = true;
            for (auto [dx, dy] : deltas)
            {
                int nx = x + dx, ny = y + dy;
                if (field(nx, ny) != '#' && last_use(nx, ny) < UT - 1 &&
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
        last_use(x, y) = UT;
        for (auto [dx, dy] : deltas)
        {
            int nx = x + dx, ny = y + dy;
            if (field(nx, ny) == '#' || last_use(nx, ny) == UT ||
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
            if (field(nx, ny) == '#' || last_use(nx, ny) == UT)
            {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < EPSILON)
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

        array<vType, deltas.size()> v;

        void swap_with(int x, int y, FluidSim* sim)
        {
            swap(sim->field(x, y), type);
            swap(sim->p(x, y), cur_p);
            swap(sim->velocity.v(x, y), v);
        }
    };

    bool propagate_move(int x, int y, bool is_first)
    {
        last_use(x, y) = UT - is_first;
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
                if (field(nx, ny) == '#' || last_use(nx, ny) == UT)
                {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < EPSILON)
                {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (abs(sum) < EPSILON)
            {
                break;
            }

            pType p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(abs(velocity.get(x, y, dx, dy)) > EPSILON &&
                   field(nx, ny) != '#' && last_use(nx, ny) < UT);

            ret = (last_use(nx, ny) == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use(x, y) = UT;
        for (size_t i = 0; i < deltas.size(); ++i)
        {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field(nx, ny) != '#' && last_use(nx, ny) < UT - 1 &&
                abs(velocity.get(x, y, dx, dy)) < EPSILON)
            {
                propagate_stop(nx, ny);
            }
        }
        if (ret)
        {
            if (!is_first)
            {
                ParticleParams pp{};
                debug_thread_matrix[x][y] = '$';
                pp.swap_with(x, y, this);
                pp.swap_with(nx, ny, this);
                pp.swap_with(x, y, this);
            }
        }
        return ret;
    }

    void save_to_last_used_buffer()
    {
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                last_use_copy(i, j) = last_use(i, j);
            }
        }
    }

    void run()
    {
        cout << "Starting simultaion..." << "\n";
#ifdef BENCH
        auto start = chrono::high_resolution_clock::now();
#endif
        rho[' '] = 0.01;
        rho['.'] = 1000;
        pType g = 0.1;

#ifdef THREAD_NUM
#pragma omp parallel for collapse(2)
#endif
        for (size_t x = 0; x < n; ++x)
        {
            for (size_t y = 0; y < m; ++y)
            {
                if (field(x, y) == '#')
                    continue;
                for (auto [dx, dy] : deltas)
                {
                    dirs(x, y) +=
                        static_cast<pType>(field(x + dx, y + dy) != '#');
                }
            }
        }

        for (size_t i = 0; i < T; ++i)
        {

            pType total_delta_p = 0;
// Apply external forces
#ifdef THREAD_NUM
#pragma omp parallel for collapse(2)
#endif
            for (size_t x = 0; x < n; ++x)
            {
                for (size_t y = 0; y < m; ++y)
                {
                    if (field(x, y) == '#')
                        continue;
                    if (field(x + 1, y) != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            // Apply forces from p
            // memcpy(sim.old_p, sim.p, sizeof(sim.p));
            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < m; ++j)
                {
                    old_p(i, j) = p(i, j);
                }
            }

            for (size_t x = 0; x < n; ++x)
            {
                for (size_t y = 0; y < m; ++y)
                {
                    if (field(x, y) == '#')
                        continue;
                    for (auto [dx, dy] : deltas)
                    {
                        int nx = x + dx, ny = y + dy;
                        if (field(nx, ny) != '#' && old_p(nx, ny) < old_p(x, y))
                        {
                            auto delta_p = old_p(x, y) - old_p(nx, ny);
                            auto force = delta_p;
                            auto& contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int)field(nx, ny)] >= force)
                            {
                                contr -= force / rho[(int)field(nx, ny)];
                                continue;
                            }
                            force -= static_cast<pType>(
                                contr * rho[(int)field(nx, ny)]);
                            contr = 0;
                            velocity.add(x, y, dx, dy,
                                         force / rho[(int)field(x, y)]);
                            p(x, y) -= force / dirs(x, y);
                            total_delta_p -= force / dirs(x, y);
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow = {};
            prop = false;
            do
            {
                UT += 2;

                prop = false;

                save_to_last_used_buffer();

                for (auto&& border : borders.borders)
                {
                    for (size_t x = 0; x < n; ++x)
                    {
                        if ((size_t)border.second < m)
                        {
                            if (field(x, border.second) != '#' &&
                                last_use(x, border.second) != UT)
                            {
                                auto [t, local_prop, _] = propagate_flow(
                                    x, border.second, 1, last_use, false);
                                if (t > EPSILON)
                                {
                                    prop = true;
                                }
                            }
                        }
                    }
                }

                start_barrier.arrive_and_wait();
                end_barrier.arrive_and_wait();

                for (size_t x = 0; x < n; ++x)
                {
                    for (size_t y = 0; y < 0; ++y)
                    {
                        last_use(x, y) =
                            max(last_use(x, y), last_use_copy(x, y));
                    }
                }

            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < n; ++x)
            {
                for (size_t y = 0; y < m; ++y)
                {
                    if (field(x, y) == '#')
                        continue;
                    for (auto [dx, dy] : deltas)
                    {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (abs(old_v) > EPSILON)
                        {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force =
                                (old_v - new_v) * rho[(int)field(x, y)];
                            if (field(x, y) == '.')
                                force *= 0.8;
                            if (field(x + dx, y + dy) == '#')
                            {
                                p(x, y) += force / dirs(x, y);
                                total_delta_p += force / dirs(x, y);
                            }
                            else
                            {
                                p(x + dx, y + dy) +=
                                    force / dirs(x + dx, y + dy);
                                total_delta_p += force / dirs(x + dx, y + dy);
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < n; ++x)
            {
                for (size_t y = 0; y < m; ++y)
                {
                    if (field(x, y) != '#' && last_use(x, y) != UT)
                    {
                        if (random01() < move_prob(x, y))
                        {
                            prop = true;
                            propagate_move(x, y, true);
                        }
                        else
                        {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            // if (prop)
            // {
            //     cout << "Tick " << i << ":\n";
            //     for (size_t x = 0; x < n; ++x)
            //     {
            //         for (size_t y = 0; y < m; ++y)
            //         {
            //             cout << debug_thread_matrix[x][y];
            //         }
            //         cout << "\n";
            //     }
            // }

            if (prop)
            {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < n; ++x)
                {
                    for (size_t y = 0; y < m; ++y)
                    {
                        cout << field(x, y);
                    }
                    cout << "\n";
                }
            }

#ifdef BENCH
            if (i >= 200)
            {
                auto now = chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = now - start;
                cout << "Time elapsed: " << elapsed.count() << " s\n";
                cout << "Waiting time elapsed " << debug_elapsed_waiting
                     << " s\n";
                cout << "Props for thread: ";
                for (auto a : debug_thread)
                {
                    cout << a << " ";
                }
                cout << "\n";
                cout << "Missing border stat: " << debug_missed.first << " "
                     << debug_missed.second << " "
                     << (double)debug_missed.first / debug_missed.second * 100
                     << "%" << "\n";
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < n; ++x)
                {
                    for (size_t y = 0; y < m; ++y)
                    {
                        cout << field(x, y);
                    }
                    cout << "\n";
                }
                return;
                // exit(0);
            }
#endif
        }
    }
};

void printDelim()
{
    cout << "#############################################" << "\n";
}

int main(int argc, char** argv)
{

    // Outpud types debug info

    printDelim();

    cout << "Types string: " << TYPES_STRING << "\n";

    cout << "Types amount: " << TYPES_COUNT << "\n";

    cout << "Type list: ";

    for (auto&& s : typesStrArr)
    {
        cout << s << " ";
    }

    cout << "\n";

    // Outpud sizes debug info

    printDelim();

    cout << "Sizes string: " << SIZES_STRING << "\n";

    cout << "Sizes amount: " << SIZES_COUNT << "\n";

    // Parse args

    string pTypeStr;
    string vTypeStr;
    string vFlowTypeStr;
    string fieldPath;
    int numThreads;
    if (argc < 4)
    {
        cout << "Not enogh arguments" << "\n";
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

    cout << "Types that will be used: " << "\n";
    cout << "p-type: " << pTypeStr << "\n";
    cout << "v-type: " << vTypeStr << "\n";
    cout << "v-flow-type: " << vFlowTypeStr << "\n";

    assert(!pTypeStr.empty());
    assert(!vTypeStr.empty());
    assert(!vFlowTypeStr.empty());

    // Create field

    printDelim();

    FieldReader fr{ fieldPath };

    global_n = fr.N;
    global_m = fr.M;

    cout << "Field read successfully" << "\n";
    cout << "Field sizes: " << fr.N << "x" << fr.M << "\n";

    // Thread num

    printDelim();

    cout << "Number of threads: " << numThreads << "\n";

    // Choose instance

    printDelim();

    bool sizeMatches = false;

    strToType(
        pTypeStr,
        [&]<typename pType>
        {
            cout << "pType chosen successfully: " << pTypeStr << "\n";
            strToType(
                vTypeStr,
                [&]<typename vType>
                {
                    cout << "vType chosen successfully: " << vTypeStr << "\n";
                    strToType(
                        vFlowTypeStr,
                        [&]<typename vFlowType>
                        {
                            cout << "vFlowType chosen successfully: "
                                 << vFlowTypeStr << "\n";
                            printDelim();
                            numsToSize(fr.N, fr.M,
                                       [&]<typename sizeType>
                                       {
                                           sizeMatches = true;
                                           cout << "Sizes matched: " << fr.N
                                                << "x" << fr.M << "\n";
                                           FluidSim<pType, vType, vFlowType,
                                                    sizeType::n, sizeType::m>
                                               sim(fr, numThreads);
                                           sim.run();
                                       });

                            if (!sizeMatches)
                            {
                                cout << "Sizes do not match" << "\n";
                                FluidSim<pType, vType, vFlowType, 0, 0> sim(
                                    fr, numThreads);
                                sim.run();
                            }
                        });
                });
        });

    return 0;
}
