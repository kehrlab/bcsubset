#include "workflow.h"

using namespace std;

int main(int argc, char const * argv[])
{
    // Check right number or arguments for execution
    if (argc < 3)
    {
        printHelp("bcsubset");
        return 1;
    }

    // Correct parsing of command line to specify commands (bcprocess and consensus)
    const char * command = argv[1];

    if (strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0)
    {
        printHelp("bcsubset");
        return 0;
    }
    else
    {
        std::cout << "\n[bcsubset]   Bam file will be processed\n" <<std::endl;
        return bamSubset(argc, argv);
    }

    return 0;
}

