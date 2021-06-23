#ifndef ARGPARSE_H_
#define ARGPARSE_H_

#include <iostream>
#include <seqan/arg_parse.h>

using namespace seqan;

struct Parameters
{
    CharString bamFileName;
    //CharString baiFileName;
    CharString bcWlFileName;
    CharString outBamFileName;

    //Parameters(): bamFileName(""), baiFileName(""), bcWlFileName(""), outBamFileName(""){}
};


// ==========================================================================
// Function printHelp()
// ==========================================================================
// Print the help lines.

void printHelp(char const * name)
{
    std::cerr << std::endl;
    std::cerr << "BCSubset - Bam file subsetting by barcode" << std::endl;
    std::cerr << "=========================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << name << " BAMFILE -f BC_WHITELIST_FILENAME -o OUT_BAMFILE_NAME " /*COMMAND\033[0m [\033[4mOPTIONS\033[0m]*/ << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mOPTIONS\033[0m" << std::endl;
    std::cerr << "    \033[1m-f\033[0m   Specify a file containing the whitelisted barcodes to use for subsetting the bam file." << std::endl;
    std::cerr << std::endl;
    std::cerr << "    \033[1m-o\033[0m   Specify an output name for the bam file containing output reads with the whitelisted barcodes." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << "BCSubset" << " version: " << VERSION << std::endl;
    std::cerr << "    Last update " << DATE << std::endl;
    std::cerr << "    Contact: Ana Pinson (ana.pinson[at]bihealth.de)" << std::endl;
    std::cerr << std::endl;
    //std::cerr << "Try `" << name << " --help' for more information on each command." << std::endl;
    //std::cerr << std::endl;
}

ArgumentParser::ParseResult parseCommandLine(Parameters & params, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("bcsubset");

    // Bam File
    addArgument(parser, ArgParseArgument(
        ArgParseArgument::STRING, "BAMFILE"));
    // BC whitelist File
    addOption(parser, ArgParseOption(
        "f", "bc_whitelist", "Location of file containing barcode whitelist",
        ArgParseArgument::INPUT_FILE, "FILE"));
    // Out Bam File Name
    addOption(parser, ArgParseOption(
        "o", "out_bamfile", "Output name for barcode subset bam file",
        ArgParseArgument::OUTPUT_FILE, "FILE"));
    
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values and print them.
    getArgumentValue(params.bamFileName, parser, 0);

    getOptionValue(params.bcWlFileName, parser, "bc_whitelist");

    getOptionValue(params.outBamFileName, parser, "out_bamfile");
    
    return ArgumentParser::PARSE_OK;
}

#endif /* ARGPARSE_H_ */