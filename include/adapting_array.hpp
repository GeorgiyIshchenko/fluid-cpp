#include <array>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

template <typename T, size_t N, size_t M, bool is_static = true>
class AdaptingArray
{
private:
    using StorageType =
        typename std::conditional<is_static, std::array<T, N * M>,
                                  std::vector<T>>::type;
    StorageType data;

    size_t n = N;
    size_t m = M;

public:
    AdaptingArray()
        requires is_static
        : data{}
    {
    }

    AdaptingArray()
        requires(!is_static)
        : data{}
    {
    }

    explicit AdaptingArray(const T& default_value)
    {
        if constexpr (is_static)
        {
            std::fill(data.begin(), data.end(), default_value);
        }
        else
        {
            data.resize(n * m, default_value);
        }
    }

    AdaptingArray(size_t rows, size_t cols, const T& default_value = T{})
        requires(!is_static): n(rows), m(cols)
    {
        data.resize(rows * cols, default_value);
    }

    T& operator()(size_t row, size_t col)
    {
        return data[row * m + col];
    }

    const T& operator()(size_t row, size_t col) const
    {
        return data[row * m + col];
    }

    constexpr size_t get_rows() const
    {
        return n;
    }
    constexpr size_t get_cols() const
    {
        return m;
    }

    void fill(const T& value)
    {
        if constexpr (is_static)
        {
            std::fill(data.begin(), data.end(), value);
        }
        else
        {
            std::fill(data.begin(), data.end(), value);
        }
    }

    void print() const
    {
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < m; ++j)
            {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
};