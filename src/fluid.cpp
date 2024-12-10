#include <array>
#include <bits/stdc++.h>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <ranges>
#include <regex>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <type_traits>

using namespace std;

// helpers

template <typename... Types> struct CountStruct;

template <typename First, typename... Others>
struct CountStruct<First, Others...>
{
    static const size_t val = 1 + CountStruct<Others...>::val;
};

template <> struct CountStruct<>
{
    static const size_t val = 0;
};

// sizes business

#ifndef DSIZES
#define DSIZES S(36, 84)
#endif

template <int N, int M> struct CompileSize
{
    static const int n = N;
    static const int m = M;
};

#define S(x, y) CompileSize<x, y>
#define SIZES_STRING STRINGIFY((SIZES))

#define SIZES_COUNT CountStruct<SIZES>::val

// types business

#define DOUBLE double
#define FLOAT float
#define FIXED(x, y) Fixed<x, y>
#define FAST_FIXED(x, y) FastFixed<x, y>

#ifndef TYPES
#define TYPES
#endif

#define STRINGIFY_IMPL(x) #x
#define STRINGIFY(x) STRINGIFY_IMPL(x)
#define TYPES_STRING STRINGIFY((TYPES))

#define TYPES_COUNT CountStruct<TYPES>::val

template <size_t N>
struct SuitableType
    : type_identity<
          conditional<N == 8, int8_t,
                      conditional<N == 16, int16_t,
                                  conditional<N == 32, int32_t,
                                              enable_if<N == 64, int64_t>>>>>
{
};

template <size_t N>
struct SuitableFastType
    : type_identity<conditional<
          N <= 8, int_fast8_t,
          conditional<N <= 16, int_fast16_t,
                      conditional<N <= 32, int_fast32_t,
                                  enable_if<N <= 64, int_fast64_t>>>>>
{
};

template <size_t N, size_t M> struct Fixed
{

    constexpr Fixed(int v)
        : v(v << M)
    {
    }
    constexpr Fixed(float f)
        : v(f * (1 << M))
    {
    }
    constexpr Fixed(double f)
        : v(f * (1 << M))
    {
    }
    constexpr Fixed()
        : v(0)
    {
    }

    static constexpr Fixed from_raw(SuitableType<N> x)
    {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    SuitableType<N> v;

    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
};

template <int N, int M> struct FastFixed
{

    constexpr FastFixed(int v)
        : v(v << M)
    {
    }
    constexpr FastFixed(float f)
        : v(f * (1 << M))
    {
    }
    constexpr FastFixed(double f)
        : v(f * (1 << M))
    {
    }
    constexpr FastFixed()
        : v(0)
    {
    }

    static constexpr FastFixed from_raw(SuitableFastType<N> x)
    {
        FastFixed ret;
        ret.v = x;
        return ret;
    }

    SuitableFastType<N> v;

    auto operator<=>(const FastFixed&) const = default;
    bool operator==(const FastFixed&) const = default;
};

constexpr array<string, TYPES_COUNT> parseTypesStr()
{
    array<string, TYPES_COUNT> result;
    string typesString(TYPES_STRING);
    size_t current_word_idx = 0;
    string curent_word = {};
    bool openedAngleBracket = false;
    for (size_t i = 1; i < typesString.size() - 1; ++i)
    {
        if (typesString[i] == '<')
        {
            openedAngleBracket = true;
        }
        else if (typesString[i] == '>')
        {
            openedAngleBracket = false;
        }

        if (typesString[i] == ',' && !openedAngleBracket)
        {
            result[current_word_idx] = curent_word;
            ++current_word_idx;
            curent_word = {};
        }
        else
        {
            curent_word += typesString[i];
        }
    }
    result[TYPES_COUNT - 1] = curent_word;
    return result;
};

auto typesStrArr = parseTypesStr();

template <typename F, typename... Ts, size_t... Is>
constexpr auto chooseTemplate(std::string_view name, F f,
                              std::index_sequence<Is...> ind) -> void
{
    (void)ind;
    (
        [&]<size_t i>
        {
            if (typesStrArr[i] == name)
            {
                f.template operator()<std::remove_cvref_t<decltype(std::get<i>(
                    std::tuple<Ts...>()))>>();
            }
        }.template operator()<Is>(),
        ...);
}

void strToType(std::string_view arg, auto f)
{
    chooseTemplate<decltype(f), TYPES>(arg, f,
                                       std::index_sequence_for<TYPES>{});
}

string strToTypeStr(string& name)
{
    if (name == "FLOAT" || name == "DOUBLE")
    {
        // convert for lower case
        std::transform(name.begin(), name.end(), name.begin(),
                       [](unsigned char c) { return std::tolower(c); });
        return name;
    }
    // fixed and fast_fixed
    const regex fixedRegex("Fixed\\((\\d+),\\s*(\\d+)\\)",
                           std::regex_constants::icase);
    const regex fastFixedRegex("Fast_Fixed\\((\\d+),\\s*(\\d+)\\)",
                               std::regex_constants::icase);
    std::smatch matches;
    if (regex_search(name, matches, fastFixedRegex))
    {
        return "FastFixed<" + matches[1].str() + "," + matches[2].str() + ">";
    }
    else if (regex_search(name, matches, fixedRegex))
    {
        return "Fixed<" + matches[1].str() + "," + matches[2].str() + ">";
    }
    throw runtime_error("Unknown input type: " + name);
}

constexpr size_t N = 36, M = 84;
// constexpr size_t N = 14, M = 5;
// constexpr size_t T = 1'000'000;
// constexpr std::array<pair<int, int>, 4> deltas{
//     { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } }
// };

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

char field[N][M + 1] = {
    "##########################################################################"
    "##########",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "#                                       ......... "
    "         #",
    "#..............#            #           ......... "
    "         #",
    "#..............#            #           ......... "
    "         #",
    "#..............#            #           ......... "
    "         #",
    "#..............#            # "
    "         #",
    "#..............#            # "
    "         #",
    "#..............#            # "
    "         #",
    "#..............#            # "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............# "
    "         #",
    "#..............#............################                     # "
    "#",
    "#...........................#....................................# "
    "#",
    "#...........................#....................................# "
    "#",
    "#...........................#....................................# "
    "#",
    "################################################################## "
    "#",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "# "
    "         #",
    "##########################################################################"
    "##########",
};

template <typename pType, typename vType, typename vFlowType, int N, int M>
class FluidSim
{
    // static constexpr Fixed inf =
    //     Fixed::from_raw(std::numeric_limits<int32_t>::max());
    // static constexpr Fixed eps = Fixed::from_raw(deltas.size());

    // Fixed operator+(Fixed a, Fixed b)
    // {
    //     return Fixed::from_raw(a.v + b.v);
    // }

    // Fixed operator-(Fixed a, Fixed b)
    // {
    //     return Fixed::from_raw(a.v - b.v);
    // }

    // Fixed operator*(Fixed a, Fixed b)
    // {
    //     return Fixed::from_raw(((int64_t)a.v * b.v) >> 16);
    // }

    // Fixed operator/(Fixed a, Fixed b)
    // {
    //     return Fixed::from_raw(((int64_t)a.v << 16) / b.v);
    // }

    // Fixed& operator+=(Fixed& a, Fixed b)
    // {
    //     return a = a + b;
    // }

    // Fixed& operator-=(Fixed& a, Fixed b)
    // {
    //     return a = a - b;
    // }

    // Fixed& operator*=(Fixed& a, Fixed b)
    // {
    //     return a = a * b;
    // }

    // Fixed& operator/=(Fixed& a, Fixed b)
    // {
    //     return a = a / b;
    // }

    // Fixed operator-(Fixed x)
    // {
    //     return Fixed::from_raw(-x.v);
    // }

    // Fixed abs(Fixed x)
    // {
    //     if (x.v < 0)
    //     {
    //         x.v = -x.v;
    //     }
    //     return x;
    // }

    // ostream& operator<<(ostream& out, Fixed x)
    // {
    //     return out << x.v / (double)(1 << 16);
    // }

    // Fixed rho[256];

    // Fixed p[N][M]{}, old_p[N][M];

    // struct VectorField
    // {
    //     array<Fixed, deltas.size()> v[N][M];
    //     Fixed& add(int x, int y, int dx, int dy, Fixed dv)
    //     {
    //         return get(x, y, dx, dy) += dv;
    //     }

    //     Fixed& get(int x, int y, int dx, int dy)
    //     {
    //         size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
    //         assert(i < deltas.size());
    //         return v[x][y][i];
    //     }
    // };

    // VectorField velocity{}, velocity_flow{};
    // int last_use[N][M]{};
    // int UT = 0;

    // mt19937 rnd(1337);

    // tuple<Fixed, bool, pair<int, int>> propagate_flow(int x, int y, Fixed lim)
    // {
    //     last_use[x][y] = UT - 1;
    //     Fixed ret = 0;
    //     for (auto [dx, dy] : deltas)
    //     {
    //         int nx = x + dx, ny = y + dy;
    //         if (field[nx][ny] != '#' && last_use[nx][ny] < UT)
    //         {
    //             auto cap = velocity.get(x, y, dx, dy);
    //             auto flow = velocity_flow.get(x, y, dx, dy);
    //             if (flow == cap)
    //             {
    //                 continue;
    //             }
    //             // assert(v >= velocity_flow.get(x, y, dx, dy));
    //             auto vp = min(lim, cap - flow);
    //             if (last_use[nx][ny] == UT - 1)
    //             {
    //                 velocity_flow.add(x, y, dx, dy, vp);
    //                 last_use[x][y] = UT;
    //                 // cerr << x << " " << y << " -> " << nx << " " << ny << " "
    //                 <<
    //                     // vp << " / " << lim << "\n";
    //                     return { vp, 1, { nx, ny } };
    //             }
    //             auto [t, prop, end] = propagate_flow(nx, ny, vp);
    //             ret += t;
    //             if (prop)
    //             {
    //                 velocity_flow.add(x, y, dx, dy, t);
    //                 last_use[x][y] = UT;
    //                 // cerr << x << " " << y << " -> " << nx << " " << ny << " "
    //                 <<
    //                     // t << " / " << lim << "\n";
    //                     return { t, prop && end != pair(x, y), end };
    //             }
    //         }
    //     }
    //     last_use[x][y] = UT;
    //     return { ret, 0, { 0, 0 } };
    // }

    // Fixed random01()
    // {
    //     return Fixed::from_raw((rnd() & ((1 << 16) - 1)));
    // }

    // void propagate_stop(int x, int y, bool force = false)
    // {
    //     if (!force)
    //     {
    //         bool stop = true;
    //         for (auto [dx, dy] : deltas)
    //         {
    //             int nx = x + dx, ny = y + dy;
    //             if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
    //                 velocity.get(x, y, dx, dy) > 0)
    //             {
    //                 stop = false;
    //                 break;
    //             }
    //         }
    //         if (!stop)
    //         {
    //             return;
    //         }
    //     }
    //     last_use[x][y] = UT;
    //     for (auto [dx, dy] : deltas)
    //     {
    //         int nx = x + dx, ny = y + dy;
    //         if (field[nx][ny] == '#' || last_use[nx][ny] == UT ||
    //             velocity.get(x, y, dx, dy) > 0)
    //         {
    //             continue;
    //         }
    //         propagate_stop(nx, ny);
    //     }
    // }

    // Fixed move_prob(int x, int y)
    // {
    //     Fixed sum = 0;
    //     for (size_t i = 0; i < deltas.size(); ++i)
    //     {
    //         auto [dx, dy] = deltas[i];
    //         int nx = x + dx, ny = y + dy;
    //         if (field[nx][ny] == '#' || last_use[nx][ny] == UT)
    //         {
    //             continue;
    //         }
    //         auto v = velocity.get(x, y, dx, dy);
    //         if (v < 0)
    //         {
    //             continue;
    //         }
    //         sum += v;
    //     }
    //     return sum;
    // }

    // struct ParticleParams
    // {
    //     char type;
    //     Fixed cur_p;
    //     array<Fixed, deltas.size()> v;

    //     void swap_with(int x, int y)
    //     {
    //         swap(field[x][y], type);
    //         swap(p[x][y], cur_p);
    //         swap(velocity.v[x][y], v);
    //     }
    // };

    // bool propagate_move(int x, int y, bool is_first)
    // {
    //     last_use[x][y] = UT - is_first;
    //     bool ret = false;
    //     int nx = -1, ny = -1;
    //     do
    //     {
    //         std::array<Fixed, deltas.size()> tres;
    //         Fixed sum = 0;
    //         for (size_t i = 0; i < deltas.size(); ++i)
    //         {
    //             auto [dx, dy] = deltas[i];
    //             int nx = x + dx, ny = y + dy;
    //             if (field[nx][ny] == '#' || last_use[nx][ny] == UT)
    //             {
    //                 tres[i] = sum;
    //                 continue;
    //             }
    //             auto v = velocity.get(x, y, dx, dy);
    //             if (v < 0)
    //             {
    //                 tres[i] = sum;
    //                 continue;
    //             }
    //             sum += v;
    //             tres[i] = sum;
    //         }

    //         if (sum == 0)
    //         {
    //             break;
    //         }

    //         Fixed p = random01() * sum;
    //         size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

    //         auto [dx, dy] = deltas[d];
    //         nx = x + dx;
    //         ny = y + dy;
    //         assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' &&
    //                last_use[nx][ny] < UT);

    //         ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    //     } while (!ret);
    //     last_use[x][y] = UT;
    //     for (size_t i = 0; i < deltas.size(); ++i)
    //     {
    //         auto [dx, dy] = deltas[i];
    //         int nx = x + dx, ny = y + dy;
    //         if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
    //             velocity.get(x, y, dx, dy) < 0)
    //         {
    //             propagate_stop(nx, ny);
    //         }
    //     }
    //     if (ret)
    //     {
    //         if (!is_first)
    //         {
    //             ParticleParams pp{};
    //             pp.swap_with(x, y);
    //             pp.swap_with(nx, ny);
    //             pp.swap_with(x, y);
    //         }
    //     }
    //     return ret;
    // }

    // int dirs[N][M]{};
};

// int main() {
//     rho[' '] = 0.01;
//     rho['.'] = 1000;
//     Fixed g = 0.1;

//     for (size_t x = 0; x < N; ++x) {
//         for (size_t y = 0; y < M; ++y) {
//             if (field[x][y] == '#')
//                 continue;
//             for (auto [dx, dy] : deltas) {
//                 dirs[x][y] += (field[x + dx][y + dy] != '#');
//             }
//         }
//     }

//     for (size_t i = 0; i < T; ++i) {

//         Fixed total_delta_p = 0;
//         // Apply external forces
//         for (size_t x = 0; x < N; ++x) {
//             for (size_t y = 0; y < M; ++y) {
//                 if (field[x][y] == '#')
//                     continue;
//                 if (field[x + 1][y] != '#')
//                     velocity.add(x, y, 1, 0, g);
//             }
//         }

//         // Apply forces from p
//         memcpy(old_p, p, sizeof(p));
//         for (size_t x = 0; x < N; ++x) {
//             for (size_t y = 0; y < M; ++y) {
//                 if (field[x][y] == '#')
//                     continue;
//                 for (auto [dx, dy] : deltas) {
//                     int nx = x + dx, ny = y + dy;
//                     if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y])
//                     {
//                         auto delta_p = old_p[x][y] - old_p[nx][ny];
//                         auto force = delta_p;
//                         auto &contr = velocity.get(nx, ny, -dx, -dy);
//                         if (contr * rho[(int) field[nx][ny]] >= force) {
//                             contr -= force / rho[(int) field[nx][ny]];
//                             continue;
//                         }
//                         force -= contr * rho[(int) field[nx][ny]];
//                         contr = 0;
//                         velocity.add(x, y, dx, dy, force / rho[(int)
//                         field[x][y]]); p[x][y] -= force / dirs[x][y];
//                         total_delta_p -= force / dirs[x][y];
//                     }
//                 }
//             }
//         }

//         // Make flow from velocities
//         velocity_flow = {};
//         bool prop = false;
//         do {
//             UT += 2;
//             prop = 0;
//             for (size_t x = 0; x < N; ++x) {
//                 for (size_t y = 0; y < M; ++y) {
//                     if (field[x][y] != '#' && last_use[x][y] != UT) {
//                         auto [t, local_prop, _] = propagate_flow(x, y, 1);
//                         if (t > 0) {
//                             prop = 1;
//                         }
//                     }
//                 }
//             }
//         } while (prop);

//         // Recalculate p with kinetic energy
//         for (size_t x = 0; x < N; ++x) {
//             for (size_t y = 0; y < M; ++y) {
//                 if (field[x][y] == '#')
//                     continue;
//                 for (auto [dx, dy] : deltas) {
//                     auto old_v = velocity.get(x, y, dx, dy);
//                     auto new_v = velocity_flow.get(x, y, dx, dy);
//                     if (old_v > 0) {
//                         assert(new_v <= old_v);
//                         velocity.get(x, y, dx, dy) = new_v;
//                         auto force = (old_v - new_v) * rho[(int)
//                         field[x][y]]; if (field[x][y] == '.')
//                             force *= 0.8;
//                         if (field[x + dx][y + dy] == '#') {
//                             p[x][y] += force / dirs[x][y];
//                             total_delta_p += force / dirs[x][y];
//                         } else {
//                             p[x + dx][y + dy] += force / dirs[x + dx][y +
//                             dy]; total_delta_p += force / dirs[x + dx][y +
//                             dy];
//                         }
//                     }
//                 }
//             }
//         }

//         UT += 2;
//         prop = false;
//         for (size_t x = 0; x < N; ++x) {
//             for (size_t y = 0; y < M; ++y) {
//                 if (field[x][y] != '#' && last_use[x][y] != UT) {
//                     if (random01() < move_prob(x, y)) {
//                         prop = true;
//                         propagate_move(x, y, true);
//                     } else {
//                         propagate_stop(x, y, true);
//                     }
//                 }
//             }
//         }

//         if (prop) {
//             cout << "Tick " << i << ":\n";
//             for (size_t x = 0; x < N; ++x) {
//                 cout << field[x] << "\n";
//             }
//         }
//     }
// }

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

    // Choose instance

    strToType(pTypeStr,
              [&]<typename pType>
              {
                  strToType(vTypeStr,
                            [&]<typename vType>
                            {
                                strToType(vFlowTypeStr,
                                          [&]<typename vFlowType>
                                          {
                                            FluidSim<pType, vType, vFlowType, N, M> sim;
                                            (void) sim;
                                          });
                            });
              });

    return 0;
}
