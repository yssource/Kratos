#pragma once

#include <vector>

namespace Kratos {

template <typename TData>
class Grid
{
public:
    using DataType = TData;

private:
    int m_size1;
    int m_size2;
    std::vector<DataType> m_values;

private:
    inline int
    Index(
        const int Index1,
        const int Index2) const noexcept
    {
        return Index1 * Size2() + Index2;
    }

public:
    Grid()
    : m_size1(0)
    , m_size2(0)
    , m_values(0)
    {
    }

    Grid(
        const int Size1,
        const int Size2)
    : m_size1(Size1)
    , m_size2(Size2)
    , m_values(Size1 * Size2)
    {
    }

    int
    Size1() const
    {
        return m_size1;
    }

    int
    Size2() const
    {
        return m_size2;
    }

    DataType&
    operator()(
        const int Index1,
        const int Index2)
    {
        int index = Index(Index1, Index2);

        return m_values[index];
    }

    const DataType&
    operator()(
        const int Index1,
        const int Index2) const
    {
        int index = Index(Index1, Index2);

        return m_values[index];
    }

    void
    Resize(
        const int Size1,
        const int Size2)
    {
        m_size1 = Size1;
        m_size2 = Size2;
        m_values.resize(Size1 * Size2);
    }
};

} // namespace ANurbs