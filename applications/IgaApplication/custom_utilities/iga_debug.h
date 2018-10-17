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
};


} // namespace Kratos

#endif // !defined(KRATOS_IGA_DEBUG_H_INCLUDED)
