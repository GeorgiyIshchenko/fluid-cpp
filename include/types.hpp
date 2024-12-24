#include <bits/stdc++.h>
#include <cstddef>

using namespace std;

// constants

#define EPSILON 1e-11

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
    using type = conditional_t<
        N == 8, int8_t,
        conditional_t<N == 16, int16_t,
                      conditional_t<N == 32, int32_t,
                                    conditional_t<N == 64, int64_t, int>>>>;
};

template <size_t N> struct SuitableFastType
{
    using type = conditional_t<
        N <= 8, int_fast8_t,
        conditional_t<
            N <= 16, int_fast16_t,
            conditional_t<N <= 32, int_fast32_t,
                          conditional_t<N <= 64, int_fast64_t, int>>>>;
};

template <size_t N, size_t M, bool isFast> struct Fixed
{
    using valueType =
        std::conditional_t<isFast, typename SuitableFastType<N>::type,
                           typename SuitableType<N>::type>;
    using typeToCast = int64_t;

    constexpr Fixed(int f)
    {
        v = f << M;
    }

    constexpr Fixed(float f)
    {
        v = f * (1ULL << M);
    }
    constexpr Fixed(double f)
    {
        v = f * (1ULL << M);
    }
    constexpr Fixed()
        : v(0)
    {
    }

    template <size_t N2, size_t M2, bool isFast2>
    constexpr Fixed(const Fixed<N2, M2, isFast2>& f)
    {
        if (M > M2)
        {
            v = f.v << (M - M2); // NOLINT
        }
        else
        {
            v = f.v >> (M2 - M); // NOLINT
        }
    }

    explicit operator double() const
    {
        return static_cast<double>(v) / (1ULL << M);
    }

    explicit operator float() const
    {
        return static_cast<float>(v) / (1ULL << M);
    }

    static constexpr Fixed from_raw(
        typename conditional_t<isFast, SuitableFastType<N>, SuitableType<N>>::type x)
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
    Fixed<N1, K1, isFast1> bToA(b);
    return a.v <=> bToA.v;
}

template <size_t N1, size_t K1, bool isFast1, size_t N2, size_t K2,
          bool isFast2>
auto operator==(const Fixed<N1, K1, isFast1>& a,
                const Fixed<N2, K2, isFast2>& b)
{
    // return fabs(static_cast<double>(a) - static_cast<double>(b)) < EPSILON;
    return a.v == b.v;
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

template <floating_point T, size_t N1, size_t K1, bool isFast1>
auto operator<=>(const Fixed<N1, K1, isFast1> a, T b)
{
    if (static_cast<T>(a) < b)
        return std::strong_ordering::less;
    if (static_cast<T>(a) > b)
        return std::strong_ordering::greater;
    return std::strong_ordering::equal;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
auto operator==(const Fixed<N1, K1, isFast1> a, T b)
{
    return fabs(static_cast<double>(a) - b) < EPSILON;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator+(Fixed<N1, K1, isFast1> a, T b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(a.v + bFixed.v);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator-(Fixed<N1, K1, isFast1> a, T b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(a.v - bFixed.v);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator*(Fixed<N1, K1, isFast1> a, T b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v * bFixed.v) >> K1);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator/(Fixed<N1, K1, isFast1> a, T b)
{
    Fixed<N1, K1, isFast1> bFixed(b);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)a.v << K1) / bFixed.v);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1>& operator+=(Fixed<N1, K1, isFast1>& a, T b)
{
    return a = a + b;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1, size_t N2,
          size_t K2, bool isFast2>
Fixed<N1, K1, isFast1>& operator-=(Fixed<N1, K1, isFast1>& a, T b)
{
    return a = a - b;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1>& operator*=(Fixed<N1, K1, isFast1>& a, T b)
{
    return a = a * b;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1>& operator/=(Fixed<N1, K1, isFast1>& a, T b)
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

template <floating_point T, size_t N1, size_t K1, bool isFast1>
auto operator==(T a, const Fixed<N1, K1, isFast1> b)
{
    return fabs(a - static_cast<T>(b)) < EPSILON;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator+(T a, Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(aFixed.v + b.v);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator-(T a, Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(aFixed.v - b.v);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator*(T a, Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)aFixed.v * b.v) >> K1);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
Fixed<N1, K1, isFast1> operator/(T a, Fixed<N1, K1, isFast1> b)
{
    Fixed<N1, K1, isFast1> aFixed(a);
    return Fixed<N1, K1, isFast1>::from_raw(((int64_t)aFixed.v << K1) / b.v);
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
floating_point auto& operator+=(T& a, Fixed<N1, K1, isFast1> b)
{
    a += static_cast<T>(b);
    return a;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
floating_point auto& operator-=(T& a, Fixed<N1, K1, isFast1> b)
{
    a -= static_cast<T>(b);
    return a;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
floating_point auto& operator*=(T& a, Fixed<N1, K1, isFast1> b)
{
    a *= static_cast<T>(b);
    return a;
}

template <floating_point T, size_t N1, size_t K1, bool isFast1>
floating_point auto& operator/=(T& a, Fixed<N1, K1, isFast1> b)
{
    a /= static_cast<T>(b);
    return a;
}

// END floating_point auto op FIXED2

constexpr array<string, TYPES_COUNT + SIZES_COUNT> parseTypesStr()
{
    array<string, TYPES_COUNT + SIZES_COUNT> result;
    string typesString(TYPES_STRING);
    typesString = typesString.substr(1, typesString.size() - 2);
    string sizesString(SIZES_STRING);
    sizesString = sizesString.substr(1, sizesString.size() - 2);
    typesString = typesString + "," + sizesString;
    size_t current_word_idx = 0;
    string curent_word = {};
    bool openedAngleBracket = false;
    for (size_t i = 0; i < typesString.size(); ++i)
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
    result[TYPES_COUNT + SIZES_COUNT - 1] = curent_word;
    return result;
};

auto typesStrArr = parseTypesStr();

template <typename F, typename... Ts, size_t... Is>
constexpr auto chooseTemplate(std::string name, F f,
                              std::index_sequence<Is...> ind,
                              size_t offset) -> void
{
    (void)ind;
    (
        [&]<size_t i>
        {
            if (typesStrArr[i + offset] == name)
            {
                f.template operator()<std::remove_cvref_t<decltype(std::get<i>(
                    std::tuple<Ts...>()))>>();
            }
        }.template operator()<Is>(),
        ...);
}

void strToType(std::string arg, auto f)
{
    chooseTemplate<decltype(f), TYPES>(arg, f, std::index_sequence_for<TYPES>{},
                                       0);
}

void numsToSize(size_t N, size_t M, auto f)
{
    auto compileSizeStr =
        "CompileSize<" + to_string(N) + ", " + to_string(M) + ">";
    chooseTemplate<decltype(f), SIZES>(
        compileSizeStr, f, std::index_sequence_for<SIZES>{}, TYPES_COUNT);
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

using PropagateFlowBorder = pair<int, int>;

struct PropagateFlowBorders{
    vector<PropagateFlowBorder> borders;

    PropagateFlowBorders(int m, int parts){
        int step = m / parts;
        for (int i = 0; i < parts + 1; ++i){
            if (i * step != min((i + 1) * step, m))
                borders.push_back({i * step, min((i + 1) * step, m)});
        }

        cout << "Borders: ";
        for(auto&& b: borders){
            cout << "{" << b.first << "," << b.second << "}" << " ";
        }
        cout << endl;
    }

};
