#ifndef ARGPARSE_H_
#define ARGPARSE_H_

#include <iostream>
#include <seqan/arg_parse.h>

using namespace seqan;

struct Parameters
{
    CharString bamFileName;
    CharString bcWlFileName;
    CharString outBamFileName;
    unsigned trimming;
    CharString bctag;
};

ArgumentParser::ParseResult parseCommandLine(Parameters & params, int argc, char const ** argv)
{
    // Setup ArgumentParser
    ArgumentParser parser("bcsubset");
    
    setShortDescription(parser, "Create a BAM file subset based on barcode whitelist");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "\\fI-w BARCODE-FILE\\fP \\fI-o OUTPUT-FILE\\fP \\fI[OPTIONS]\\fP \\fIBAM-FILE\\fP");

    addDescription(parser, "Selects records from the BAM file that match the barcodes provided in a whitelist.");

    // Input BAM file
    addArgument(parser, ArgParseArgument(
        ArgParseArgument::STRING, "BAMFILE"));
    // Whitelisted barcodes ile
    addOption(parser, ArgParseOption(
        "w", "whitelist", "File containing whitelisted barcodes. One barcode per line.",
        ArgParseArgument::INPUT_FILE, "FILE"));
    setRequired(parser, "w");
    // Out BAM file name
    addOption(parser, ArgParseOption(
        "o", "out", "Output name for barcode subset BAM file.",
        ArgParseArgument::OUTPUT_FILE, "FILE"));
    setRequired(parser, "o");
    // Trimming barcode
    addOption(parser, ArgParseOption(
        "t", "trim_suffix", "Trim the last n characters from barcode in input BAM file.",
        ArgParseArgument::INTEGER, "NUM"));
    addDefaultValue(parser, "t", 0);
    //setMinValue(parser, "t", 0); // line that breaks the help
    // Specify tag for barcode
    addOption(parser, ArgParseOption(
        "b", "barcode_tag", "BAM record tag containing the barcodes to compare with the whitelist.",
        ArgParseArgument::STRING, "TAG"));
    addDefaultValue(parser, "b", "CB");
    
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values
    getArgumentValue(params.bamFileName, parser, 0);

    getOptionValue(params.bcWlFileName, parser, "whitelist");

    getOptionValue(params.outBamFileName, parser, "out");

    getOptionValue(params.trimming, parser, "trim_suffix");

    getOptionValue(params.bctag, parser, "barcode_tag");
    
    return ArgumentParser::PARSE_OK;
}

inline int checkParser(const ArgumentParser::ParseResult & res)
{
    if (res == ArgumentParser::PARSE_HELP ||
        res == ArgumentParser::PARSE_VERSION ||
        res == ArgumentParser::PARSE_WRITE_CTD ||
        res == ArgumentParser::PARSE_EXPORT_HELP)
        return 0;
    else if (res != ArgumentParser::PARSE_OK)
        return 1;
    else
        return -1;
}

#endif /* ARGPARSE_H_ */