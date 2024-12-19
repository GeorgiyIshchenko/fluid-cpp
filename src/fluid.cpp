#include <array>
#include <bits/stdc++.h>
#include <cassert>
#include <cctype>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string_view>
#include <type_traits>

using namespace std;

// constants

#define EPSILON 1e-8

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
#define FIXED(x, y) Fixed<x, y, false>
#define FAST_FIXED(x, y) Fixed<x, y, true>

#ifndef TYPES
#define TYPES
#endif

#define STRINGIFY_IMPL(x) #x
#define STRINGIFY(x) STRINGIFY_IMPL(x)
#define TYPES_STRING STRINGIFY((TYPES))

#define TYPES_COUNT CountStruct<TYPES>::val

template <size_t N> struct SuitableType
{
    using type = conditional<
        N == 8, int8_t,
        conditional_t<N == 16, int16_t,
                      conditional_t<N == 32, int32_t,
                                    conditional_t<N == 64, int64_t, int>>>>::
        type;
};

template <size_t N> struct SuitableFastType
{
    using type = conditional<
        N <= 8, int_fast8_t,
        conditional_t<
            N <= 16, int_fast16_t,
            conditional_t<N <= 32, int_fast32_t,
                          conditional_t<N <= 64, int_fast64_t, int>>>>::type;
};

template <size_t N, size_t M, bool isFast> struct Fixed
{
    using valueType =
        std::conditional_t<isFast, typename SuitableFastType<N>::type,
                           typename SuitableType<N>::type>;
    using typeToCast = int64_t;

    constexpr Fixed(int f)
    {
        v = f * (1 << M);
    }

    constexpr Fixed(float f)
    {
        v = f * (1 << M);
    }
    constexpr Fixed(double f)
    {
        v = f * (1 << M);
    }
    constexpr Fixed()
        : v(0)
    {
    }

    template <size_t N2, size_t M2, bool isFast2>
    constexpr Fixed(const Fixed<N2, M2, isFast2>& f)
    {
        if (M >= M2)
        {
            v = f.v << (M - M2);
        }
        else {
            v = f.v >> (M2 - M);
        }
    }

    explicit operator double() const
    {
        return static_cast<double>(v) / (1 << M);
    }

    explicit operator float() const
    {
        return static_cast<float>(v) / (1 << M);
    }

    static constexpr Fixed from_raw(
        conditional_t<isFast, SuitableFastType<N>, SuitableType<N>>::type x)
    {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    valueType v;

    // bool operator==(const Fixed&) const = default;
};

// FIXED to FIXED

// ???

// FIXED1 op FIXED2

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
auto operator<=>(const Fixed<N1, K1, isFast1>& a,
                 const Fixed<N2, K2, isFast2>& b)
{
    if (static_cast<double>(a) < static_cast<double>(b))
        return std::strong_ordering::less;
    if (static_cast<double>(a) > static_cast<double>(b))
        return std::strong_ordering::greater;
    return std::strong_ordering::equal;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
auto operator==(const Fixed<N1, K1, isFast1>& a,
                const Fixed<N2, K2, isFast2>& b)
{
    return fabs(static_cast<double>(a) - static_cast<double>(b)) < EPSILON;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1> operator+(Fixed<N1, K1, isFast1> a,
                                 Fixed<N2, K2, isFast2> b)
{
    Fixed<N1, K1, isFast1> bToA(b);
    return Fixed<N1, K1, isFast1>::from_raw(a.v + bToA.v);
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1> operator-(Fixed<N1, K1, isFast1> a,
                                 Fixed<N2, K2, isFast2> b)
{
    Fixed<N1, K1, isFast1> bToA(b);
    return Fixed<N1, K1, isFast1>::from_raw(a.v - bToA.v);
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1> operator*(Fixed<N1, K1, isFast1> a,
                                 Fixed<N2, K2, isFast2> b)
{
    Fixed<N1, K1, isFast1> bToA(b);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v * bToA.v) >> K1);
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1> operator/(Fixed<N1, K1, isFast1> a,
                                 Fixed<N2, K2, isFast2> b)
{
    Fixed<N1, K1, isFast1> bToA(b);
    if (bToA.v == 0){
        return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v << K1) / 1);
    }
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v << K1) / bToA.v);
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1>& operator+=(Fixed<N1, K1, isFast1>& a,
                                   Fixed<N2, K2, isFast2> b)
{
    return a = a + b;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1>& operator-=(Fixed<N1, K1, isFast1>& a,
                                   Fixed<N2, K2, isFast2> b)
{
    return a = a - b;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1>& operator*=(Fixed<N1, K1, isFast1>& a,
                                   Fixed<N2, K2, isFast2> b)
{
    return a = a * b;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1>& operator/=(Fixed<N1, K1, isFast1>& a,
                                   Fixed<N2, K2, isFast2> b)
{
    return a = a / b;
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator-(Fixed<N1, K1, isFast1> x)
{
    return Fixed<N1, K1, isFast1>::from_raw(-x.v);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> abs(Fixed<N1, K1, isFast1> x)
{
    if (x.v < 0)
    {
        x.v = -x.v;
    }
    return x;
}

template <size_t N1, size_t K1, bool isFast1>
ostream& operator<<(ostream& out, Fixed<N1, K1, isFast1> x)
{
    return out << x.v / (double)(1 << K1);
}

// END FIXED1 op FIXED2

// FIXED1 op floating point

template <size_t N1, size_t K1, bool isFast1>
auto operator<=>(const Fixed<N1, K1, isFast1> a, floating_point auto b)
{
    if (static_cast<double>(a) < b)
        return std::strong_ordering::less;
    if (static_cast<double>(a) > b)
        return std::strong_ordering::greater;
    return std::strong_ordering::equal;
}

template <size_t N1, size_t K1, bool isFast1>
auto operator==(const Fixed<N1, K1, isFast1> a, floating_point auto b)
{
    return fabs(static_cast<double>(a) - b) < EPSILON;
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator+(Fixed<N1, K1, isFast1> a,
                                 floating_point auto b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(a.v + bFixed.v);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator-(Fixed<N1, K1, isFast1> a,
                                 floating_point auto b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(a.v - bFixed.v);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator*(Fixed<N1, K1, isFast1> a,
                                 floating_point auto b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v * bFixed.v) >> K1);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator/(Fixed<N1, K1, isFast1> a,
                                 floating_point auto b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v << K1) / bFixed.v);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1>& operator+=(Fixed<N1, K1, isFast1>& a,
                                   floating_point auto b)
{
    return a = a + b;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
Fixed<N1, K1, isFast1>& operator-=(Fixed<N1, K1, isFast1>& a,
                                   floating_point auto b)
{
    return a = a - b;
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1>& operator*=(Fixed<N1, K1, isFast1>& a,
                                   floating_point auto b)
{
    return a = a * b;
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1>& operator/=(Fixed<N1, K1, isFast1>& a,
                                   floating_point auto b)
{
    return a = a / b;
}

// END FIXED1 op floating_point auto

// floating point op FIXED1

template <size_t N1, size_t K1, bool isFast1>
auto operator<=>(floating_point auto a, const Fixed<N1, K1, isFast1> b)
{
    if (a < static_cast<double>(b))
        return std::strong_ordering::less;
    if (a > static_cast<double>(b))
        return std::strong_ordering::greater;
    return std::strong_ordering::equal;
}

template <size_t N1, size_t K1, bool isFast1>
auto operator==(floating_point auto a, const Fixed<N1, K1, isFast1> b)
{
    return fabs(a - static_cast<double>(b)) < EPSILON;
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator+(floating_point auto a,
                                 Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(aFixed.v + b.v);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator-(floating_point auto a,
                                 Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(aFixed.v - b.v);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator*(floating_point auto a,
                                 Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)aFixed.v * b.v) >> K1);
}

template <size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator/(floating_point auto a,
                                 Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)aFixed.v << K1) / b.v);
}

template <size_t N1, size_t K1, bool isFast1>
floating_point auto& operator+=(floating_point auto& a,
                                Fixed<N1, K1, isFast1> b)
{
    a += static_cast<double>(b);
    return a;
}

template <size_t N1, size_t K1, bool isFast1>
floating_point auto& operator-=(floating_point auto& a,
                                Fixed<N1, K1, isFast1> b)
{
    a -= static_cast<double>(b);
    return a;
}

template <size_t N1, size_t K1, bool isFast1>
floating_point auto& operator*=(floating_point auto& a,
                                Fixed<N1, K1, isFast1> b)
{
    a *= static_cast<double>(b);
    return a;
}

template <size_t N1, size_t K1, bool isFast1>
floating_point auto& operator/=(floating_point auto& a,
                                Fixed<N1, K1, isFast1> b)
{
    a /= static_cast<double>(b);
    return a;
}

// END floating_point auto op FIXED2

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
        return "Fixed<" + matches[1].str() + ", " + matches[2].str() +
               ", true>";
    }
    else if (regex_search(name, matches, fixedRegex))
    {
        return "Fixed<" + matches[1].str() + ", " + matches[2].str() +
               ", false>";
    }
    throw runtime_error("Unknown input type: " + name);
}

template <typename pType, typename vType, typename vFlowType, int N, int M>
class FluidSim
{
public:
    // constexpr size_t N = 14, M = 5;
    constexpr static size_t T = 1'000'000;
    inline static std::array<pair<int, int>, 4> deltas = {
        { { -1, 0 }, { 1, 0 }, { 0, -1 }, { 0, 1 } }
    };

    // inline static char field[N][M + 1] = {
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

    inline static char field[N][M + 1] = {
        "######################################################################"
        "##############",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                       .........                     "
        "             #",
        "#..............#            #           .........                     "
        "             #",
        "#..............#            #           .........                     "
        "             #",
        "#..............#            #           .........                     "
        "             #",
        "#..............#            #                                         "
        "             #",
        "#..............#            #                                         "
        "             #",
        "#..............#            #                                         "
        "             #",
        "#..............#            #                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............#                                         "
        "             #",
        "#..............#............################                     #    "
        "             #",
        "#...........................#....................................#    "
        "             #",
        "#...........................#....................................#    "
        "             #",
        "#...........................#....................................#    "
        "             #",
        "##################################################################    "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "#                                                                     "
        "             #",
        "######################################################################"
        "##############",
    };

    inline static vType rho[256]{};

    inline static pType p[N][M]{}, old_p[N][M]{};

    template <typename T> struct VectorField
    {
        array<T, deltas.size()> v[N][M];

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
    inline static int last_use[N][M]{};
    inline static int UT = 0;

    mt19937 rnd{ 1337 };

    FluidSim()
    {
    }

    tuple<vFlowType, bool, pair<int, int>> propagate_flow(int x, int y,
                                                          vFlowType lim)
    {
        last_use[x][y] = UT - 1;
        vFlowType ret = 0;
        for (auto [dx, dy] : deltas)
        {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT)
            {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap)
                {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(static_cast<double>(lim),
                              static_cast<double>(cap - flow));
                if (last_use[nx][ny] == UT - 1)
                {
                    velocity_flow.add(x, y, dx, dy, static_cast<double>(vp));
                    // last_use[x][y] = UT;
                    // cerr << "A " << x << " " << y << " -> " << nx << " " << ny
                    //      << " " << vp << " / " << lim << "\n";
                    return { static_cast<vFlowType>(vp), 1, { nx, ny } };
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop)
                {
                    velocity_flow.add(x, y, dx, dy, t);
                    // last_use[x][y] = UT;
                    // cerr << "B " << x << " " << y << " -> " << nx << " " << ny
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
        array<vType, deltas.size()> v;

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

    pType dirs[N][M]{};
};

template <typename pType, typename vType, typename vFlowType, size_t N,
          size_t M>
void do_main(FluidSim<pType, vType, vFlowType, N, M>& sim)
{
    sim.rho[' '] = 0.01;
    sim.rho['.'] = 1000;
    pType g = 0.1;

    for (size_t x = 0; x < N; ++x)
    {
        for (size_t y = 0; y < M; ++y)
        {
            if (sim.field[x][y] == '#')
                continue;
            for (auto [dx, dy] : sim.deltas)
            {
                sim.dirs[x][y] +=
                    static_cast<double>(sim.field[x + dx][y + dy] != '#');
            }
        }
    }

    for (size_t i = 0; i < sim.T; ++i)
    {

        pType total_delta_p = 0;
        // Apply external forces
        for (size_t x = 0; x < N; ++x)
        {
            for (size_t y = 0; y < M; ++y)
            {
                if (sim.field[x][y] == '#')
                    continue;
                if (sim.field[x + 1][y] != '#')
                    sim.velocity.add(x, y, 1, 0, g);
            }
        }

        // Apply forces from p
        memcpy(sim.old_p, sim.p, sizeof(sim.p));
        for (size_t x = 0; x < N; ++x)
        {
            for (size_t y = 0; y < M; ++y)
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
                        force -= static_cast<pType>(static_cast<double>(
                            contr * sim.rho[(int)sim.field[nx][ny]]));
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
        bool prop = false;
        do
        {
            sim.UT += 2;
            prop = 0;
            for (size_t x = 0; x < N; ++x)
            {
                for (size_t y = 0; y < M; ++y)
                {
                    if (sim.field[x][y] != '#' && sim.last_use[x][y] != sim.UT)
                    {
                        auto [t, local_prop, _] = sim.propagate_flow(x, y, 1);
                        if (t > 0.)
                        {
                            prop = 1;
                        }
                    }
                }
            }
        } while (prop);

        // Recalculate p with kinetic energy
        for (size_t x = 0; x < N; ++x)
        {
            for (size_t y = 0; y < M; ++y)
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
        prop = false;
        for (size_t x = 0; x < N; ++x)
        {
            for (size_t y = 0; y < M; ++y)
            {
                if (sim.field[x][y] != '#' && sim.last_use[x][y] != sim.UT)
                {
                    if (sim.random01() < sim.move_prob(x, y))
                    {
                        prop = true;
                        sim.propagate_move(x, y, true);
                    }
                    else
                    {
                        sim.propagate_stop(x, y, true);
                    }
                }
            }
        }

        if (prop)
        {
            cout << "Tick " << i << ":\n";
            for (size_t x = 0; x < N; ++x)
            {
                cout << sim.field[x] << "\n";
            }
        }
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

    const int N = 36;
    const int M = 84;

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
                    strToType(vFlowTypeStr,
                              [&]<typename vFlowType>
                              {
                                  cout << "vFlowType chosen successfully: "
                                       << vFlowTypeStr << endl;
                                  printDelim();
                                  cout << "Starting simultaion..." << endl;
                                  FluidSim<pType, vType, vFlowType, N, M> sim;
                                  do_main<pType, vType, vFlowType, N, M>(sim);
                              });
                });
        });

    return 0;
}
