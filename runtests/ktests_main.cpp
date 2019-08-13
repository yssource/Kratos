#include "includes/kernel.h"
#include "testing/tester.h"

int main(int argc, char *argv[])
{
    Kratos::Kernel mKernel;
    mKernel.Initialize();

    int error_code = Kratos::Testing::Tester::RunAllTestCases();

    if(error_code != 0)
    {
        std::cout << "KRATOS TERMINATED WITH ERROR" << std::endl;
        return 1;
    }
    else
    {
        std::cout << "KRATOS TERMINATED CORRECTLY" << std::endl;
        return 0;
    }
}
