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

    //Parameters(): bamFileName(""), baiFileName(""), bcWlFileName(""), outBamFileName(""){}
};

ArgumentParser::ParseResult parseCommandLine(Parameters & params, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("bcsubset");
    
    setShortDescription(parser, "Create a BAM file subset based on barcode whitelist");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addUsageLine(parser, "\\fI-w BARCODE-FILE\\fP \\fI-o OUTPUT-FILE\\fP \\fI[OPTIONS]\\fP \\fIBAM-FILE\\fP");

    addDescription(parser, "Selects records from the bam file that match the barcodes provided in a whitelist.");

    // Bam File
    addArgument(parser, ArgParseArgument(
        ArgParseArgument::STRING, "BAMFILE"));
    // BC whitelist File
    addOption(parser, ArgParseOption(
        "w", "whitelist", "File containing whitelisted barcodes. One barcode per line.",
        ArgParseArgument::INPUT_FILE, "FILE"));
    setRequired(parser, "w");
    // Out Bam File Name
    addOption(parser, ArgParseOption(
        "o", "out", "Output name for barcode subset bam file.",
        ArgParseArgument::OUTPUT_FILE, "FILE"));
    setRequired(parser, "o");
    // Trimming option
    addOption(parser, ArgParseOption(
        "t", "trim_suffix", "Trim the last n characters from barcode in bam file.",
        ArgParseArgument::INTEGER, "NUM"));
    addDefaultValue(parser, "t", 0);
    //setMinValue(parser, "t", 0); // line that breaks the help
    // Specify tag for barcode
    addOption(parser, ArgParseOption(
        "b", "barcode_tag", "Bam record tag containing the barcodes to compare with the whitelist.",
        ArgParseArgument::STRING, "TAG"));
    addDefaultValue(parser, "b", "CB");
    
    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    
    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values and print them.
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