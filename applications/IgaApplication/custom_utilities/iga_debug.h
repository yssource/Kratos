/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if !defined(KRATOS_IGA_DEBUG_H_INCLUDED)
#define KRATOS_IGA_DEBUG_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// External includes

// Project includes


namespace Kratos
{

struct IgaDebug
{
    template <typename TScalar>
    static bool CheckWithAbsoluteTolerance(
        const TScalar A,
        const TScalar B,
        const TScalar Tolerance)
    {
        if (Tolerance < 0) {
            return true;
        }

        return std::abs(A - B) < Tolerance;
    }

    static void CheckDouble(
        Parameters& ExpectedData,
        const std::string& key,
        const double ActualValue)
    {
        const auto expected_value = ExpectedData[key].GetDouble();

        if (!CheckWithAbsoluteTolerance(expected_value, ActualValue, 1e-7)) {
            std::cout << "IgaDebug: Value mismatch '" << key << "'\n"
                      << "  expected = " << expected_value << "\n"
                      << "  actual   = " << ActualValue << std::endl;

            throw std::runtime_error("IgaDebug: Value mismatch");
        }
    }

    template <typename TVector>
    static void CheckVector(
        Parameters& ExpectedData,
        const std::string& key,
        const TVector& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetVector();

        if (expected_value.size() != ActualValue.size()) {
            std::cout << "IgaDebug: " << key << ".size()\n"
                      << "  expected = " << expected_value.size() << "\n"
                      << "  actual   = " << ActualValue.size() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        for (std::size_t i = 0; i < expected_value.size(); i++) {
            if (!CheckWithAbsoluteTolerance(expected_value[i], ActualValue[i], 1e-7)) {
                std::cout << "IgaDebug: " << key << "[" << i << "]\n"
                          << "  expected = " << expected_value[i] << "\n"
                          << "  actual   = " << ActualValue[i] << std::endl;

                throw std::runtime_error("IgaDebug: Value mismatch");
            }
        }
    }

    template <typename TMatrix>
    static void CheckMatrix(
        Parameters& ExpectedData,
        const std::string& key,
        const TMatrix& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetMatrix();

        if (expected_value.size1() != ActualValue.size1()) {
            std::cout << "IgaDebug: " << key << ".size1()\n"
                      << "  expected = " << expected_value.size1() << "\n"
                      << "  actual   = " << ActualValue.size1() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        if (expected_value.size2() != ActualValue.size2()) {
            std::cout << "IgaDebug: " << key << ".size2()\n"
                      << "  expected = " << expected_value.size2() << "\n"
                      << "  actual   = " << ActualValue.size2() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        for (std::size_t i = 0; i < expected_value.size1(); i++) {
            for (std::size_t j = 0; j < expected_value.size2(); j++) {
                if (!CheckWithAbsoluteTolerance(expected_value(i, j), ActualValue(i, j), 1e-7)) {
                    std::cout << "IgaDebug: " << key << "(" << i << ", " << j << ")\n"
                              << "  expected = " << expected_value(i, j) << "\n"
                              << "  actual   = " << ActualValue(i, j) << std::endl;

                    throw std::runtime_error("IgaDebug: Value mismatch");
                }
            }
        }
    }

    template <typename TMatrix>
    static void CheckLowerMatrix(
        Parameters& ExpectedData,
        const std::string& key,
        const TMatrix& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetMatrix();

        if (expected_value.size1() != ActualValue.size1()) {
            std::cout << "IgaDebug: " << key << ".size1()\n"
                      << "  expected = " << expected_value.size1() << "\n"
                      << "  actual   = " << ActualValue.size1() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        if (expected_value.size2() != ActualValue.size2()) {
            std::cout << "IgaDebug: " << key << ".size2()\n"
                      << "  expected = " << expected_value.size2() << "\n"
                      << "  actual   = " << ActualValue.size2() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        for (std::size_t i = 0; i < expected_value.size1(); i++) {
            for (std::size_t j = 0; j <= i; j++) {
                if (!CheckWithAbsoluteTolerance(expected_value(i, j), ActualValue(i, j), 1e-7)) {
                    std::cout << "IgaDebug: " << key << "(" << i << ", " << j << ")\n"
                              << "  expected = " << expected_value(i, j) << "\n"
                              << "  actual   = " << ActualValue(i, j) << std::endl;

                    throw std::runtime_error("IgaDebug: Value mismatch");
                }
            }
        }
    }

    template <typename TGrid>
    static void CheckLowerGridComponent(
        Parameters& ExpectedData,
        const std::string& key,
        const TGrid& ActualValue,
        const int index
    )
    {
        const auto expected_value = ExpectedData[key].GetMatrix();

        if (expected_value.size1() != ActualValue.Size1()) {
            std::cout << "IgaDebug: " << key << ".size1()\n"
                      << "  expected = " << expected_value.size1() << "\n"
                      << "  actual   = " << ActualValue.Size1() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        if (expected_value.size2() != ActualValue.Size2()) {
            std::cout << "IgaDebug: " << key << ".size2()\n"
                      << "  expected = " << expected_value.size2() << "\n"
                      << "  actual   = " << ActualValue.Size2() << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        for (std::size_t i = 0; i < expected_value.size1(); i++) {
            for (std::size_t j = 0; j <= i; j++) {
                if (!CheckWithAbsoluteTolerance(expected_value(i, j), ActualValue(i, j)[index], 1e-7)) {
                    std::cout << "IgaDebug: " << key << "(" << i << ", " << j << ")[" << index << "]\n"
                              << "  expected = " << expected_value(i, j) << "\n"
                              << "  actual   = " << ActualValue(i, j)[index] << std::endl;

                    throw std::runtime_error("IgaDebug: Value mismatch");
                }
            }
        }
    }

    template <typename TVector>
    static void CheckVectorVar(
        Parameters& ExpectedData,
        const std::string& key,
        const std::vector<TVector>& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetVector();

        if (expected_value.size() != ActualValue.size() * 3) {
            std::cout << "IgaDebug: " << key << ".size()\n"
                      << "  expected = " << expected_value.size() << "\n"
                      << "  actual   = " << ActualValue.size() * 3 << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        for (std::size_t i = 0; i < ActualValue.size(); i++) {
            for (std::size_t j = 0; j < 3; j++) {
                if (!CheckWithAbsoluteTolerance(expected_value(j * ActualValue.size() + i), ActualValue[i][j], 1e-7)) {
                    std::cout << "IgaDebug: " << key << "(" << i << ", " << j << ")\n"
                              << "  expected = " << expected_value(j * ActualValue.size() + i) << "\n"
                              << "  actual   = " << ActualValue[i][j] << std::endl;

                    throw std::runtime_error("IgaDebug: Value mismatch");
                }
            }
        }
    }

    template <typename TMatrix>
    static void CheckMatrixVar(
        Parameters& ExpectedData,
        const std::string& key,
        const std::vector<TMatrix>& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetMatrix();

        int ndofs = ActualValue.size();

        for (std::size_t r = 0; r < ndofs; r++) {
        for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
            int a;
            int b;
            MatrixIndex3To2(ndofs, r, i, j, a, b);

            auto exp = expected_value(a, b);
            auto act = ActualValue[r](i, j);

            if (!CheckWithAbsoluteTolerance(exp, act, 1e-7)) {
                std::cout << "IgaDebug: " << key << "(" << r << ", " << i <<  ", " << j << ")\n"
                        << "  expected = " << exp << "\n"
                        << "  actual   = " << act << std::endl;

                throw std::runtime_error("IgaDebug: Value mismatch");
            }
        }
        }
        }
    }

    template <typename TMatrix>
    static void CheckMatrixVarVar(
        Parameters& ExpectedData,
        const std::string& key,
        const std::vector<TMatrix>& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetMatrix();

        int ndofs = sqrt(ActualValue.size());

        for (std::size_t r = 0; r < ndofs; r++) {
        for (std::size_t s = 0; s < ndofs; s++) {
        for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
            int a;
            int b;
            MatrixIndex4To2(ndofs, r, s, i, j, a, b);

            auto exp = expected_value(a, b);
            auto act = ActualValue[r * ndofs + s](i, j);

            if (!CheckWithAbsoluteTolerance(exp, act, 1e-7)) {
                std::cout << "IgaDebug: " << key << "(" << r << ", " << i <<  ", " << j << ")\n"
                        << "  expected = " << exp << "\n"
                        << "  actual   = " << act << std::endl;

                throw std::runtime_error("IgaDebug: Value mismatch");
            }
        }
        }
        }
        }
    }

    template <typename TVector>
    static void CheckVectorVarVar(
        Parameters& ExpectedData,
        const std::string& key,
        const std::vector<TVector>& ActualValue
    )
    {
        const auto expected_value = ExpectedData[key].GetMatrix();

        if (expected_value.size1() * expected_value.size2() != ActualValue.size() * 3) {
            std::cout << "IgaDebug: " << key << ".size()\n"
                      << "  expected = " << expected_value.size1() * expected_value.size2() << "\n"
                      << "  actual   = " << ActualValue.size() * 3 << std::endl;

            throw std::runtime_error("IgaDebug: Size mismatch");
        }

        int ndofs = sqrt(ActualValue.size());

        for (std::size_t r = 0; r < ndofs; r++) {
            for (std::size_t s = 0; s < ndofs; s++) {
                for (std::size_t t = 0; t < 3; t++) {
                    int u;
                    int v;
                    Index3To2(ndofs, r, s, t, u, v);

                    auto exp = expected_value(u, v);
                    auto act = ActualValue[r * ndofs + s][t];

                    if (!CheckWithAbsoluteTolerance(exp, act, 1e-7)) {
                        std::cout << "IgaDebug: " << key << "(" << r << ", " << s <<  ", " << t << ")\n"
                                << "  expected = " << exp << "\n"
                                << "  actual   = " << act << std::endl;

                        throw std::runtime_error("IgaDebug: Value mismatch");
                    }
                }
            }
        }
    }

    static void MatrixIndex3To2(
        int ndofs,
        int r, int i, int j,
        int& a, int& b)
    {
        a = i * ndofs + r;
        b = j;
    }

    static void MatrixIndex4To2(
        int ndofs,
        int r, int s, int i, int j,
        int& a, int& b)
    {
        a = i * ndofs + r;
        b = j * ndofs + s;
    }

    static void MatrixIndex2To3(
        int ndofs,
        int a, int b,
        int& r, int& t, int& u)
    {
        r = a % ndofs;
        t = a / ndofs;
        u = b;
    }

    static void Index3To2(int ndofs, int r, int s, int t, int& u, int& v)
    {
        u = t * ndofs + r;
        v = s;
    }

    static void Index2To3(int ndofs, int u, int v, int& r, int& s, int& t)
    {
        r = u % ndofs;
        s = v;
        t = u / ndofs;
    }
};


} // namespace Kratos

#endif // !defined(KRATOS_IGA_DEBUG_H_INCLUDED)
