#ifndef BARCODE_H_
#define BARCODE_H_

#include<seqan/sequence.h>
#include<seqan/bam_io.h>
 
using namespace seqan;

// Process bam header, add @PG line
inline void processHeader(BamHeader & header, BamFileOut & bamFileOut, const CharString & command, char const ** argv)
{

  // Modify @PG
  BamHeaderRecord hrecord;
  hrecord.type = BamHeaderRecordType::BAM_HEADER_PROGRAM;
  
  // Set tags
  appendValue(hrecord.tags, Pair<CharString>("ID", "bcsubset"));
  appendValue(hrecord.tags, Pair<CharString>("VN", "0.0.1"));

  CharString clstring;
  append(clstring, "bcsubset ");
  append(clstring, command);

  for (unsigned j=1; j<length(argv);++j)
  {
    appendValue(clstring, ' ');
    append(clstring, argv[j]);
  }

  //appendValue(hrecord.tags, Pair<CharString>("CL", "ctprocess bcprocess test"));
  appendValue(hrecord.tags, Pair<CharString>("CL", clstring));
  
  appendValue(header, hrecord);

  writeHeader(bamFileOut, header);
}

/*
inline void processORHeader(BamHeader & header, BamFileOut & outReadsBamFileOut)
{
  // Modify @CO
  BamHeaderRecord hcrecord; //header comment record
  hcrecord.type = BamHeaderRecordType::BAM_HEADER_COMMENT;
  //Add encoding of reasons for reads that are filtered out in consensus
  appendValue(hcrecord.tags, Pair<CharString>("OR:", "Reasons for discarded reads OR, 0: Good reverse read, no good forward read found; 1: Not a good reverse read, nor a good forward; 2: Good reverse read in region+10K, no good forward read found in region; 3: Good forward read, no good reverse read found."));
  appendValue(header, hcrecord);

  writeHeader(outReadsBamFileOut, header);
}
*/
/*
// Extract barcode from qname in a bam record

inline void getBarcode(CharString & barcode, BamAlignmentRecord & record, unsigned barcodeLength)
{
  //std::cout<< "Barcode Length: " << barcodeLength <<std::endl;
  //qnsize string, take length
  unsigned qnlen = length(record.qName);
  CharString qname = record.qName;
  toUpper(qname);
  */

 //For bcsubset, barcodes are already tags, so only need to access tags
/*
  //Substract 5 from qnsize, and move from position qnsize-5 to end of string, barcode size can be different, not always 5 
  for (unsigned i = (qnlen - barcodeLength); i < qnlen; ++i)
  {
    if(qname[i] == 'A' || qname[i] == 'T' || qname[i] == 'C' || qname[i] == 'G')
      appendValue(barcode, qname[i]);
  }

  //std::cout<< "Barcode: " << barcode << std::endl;

  resize(record.qName, qnlen - barcodeLength -1);
}
*/

inline void getBarcodeFromTags(CharString & barcode, const BamAlignmentRecord & record)
{
  unsigned idx = 0;

  BamTagsDict tagsDict(record.tags);
  bool keyFound = findTagKey(idx, tagsDict, "CR");
  if (keyFound)
  {
    if (!extractTagValue(barcode, tagsDict, idx))
    {
      std::cerr << "ERROR: There was an error extracting BC from tags of record: " << record.qName << "\n";
    }
  }
  else
  {
    std::cerr << "ERROR: BC not found in tags of records: " << record.qName << "\n";
  }
  
}

// Add extracted barcode as a tag to all records
 
void processBarcodes(const CharString & bamFileName, const CharString & outBamFileName, char const ** argv, const unsigned barcodeLength = 5)
//void processBarcodes(const CharString & bamFileName, const CharString & outBamFileName, char const ** argv, const unsigned barcodeLength = 15)
{
  // Open input file, BamFileIn can read SAM and BAM files.
  BamFileIn bamFileIn(toCString(bamFileName));
 
  // Open output file, BamFileOut accepts also an ostream and a format tag.
  std::fstream outStream;
  outStream.open(toCString(outBamFileName), std::ios::out);
  if (!outStream.good())
    SEQAN_THROW(FileOpenError(toCString(outBamFileName)));
    
  BamFileOut bamFileOut(context(bamFileIn), outStream, Bam());
  open(bamFileOut, outStream, Bam());

  BamHeader header;
  readHeader(header, bamFileIn);
 
  // Access header
  processHeader(header, bamFileOut, "bcprocess", argv);
  
  // Access records
  BamAlignmentRecord record;
  while (!atEnd(bamFileIn))
  {
    readRecord(record, bamFileIn);

    CharString barcode;
    reserve(barcode, barcodeLength, Exact());
    getBarcode(barcode, record, barcodeLength);
 
    // Access records tags, add tag BC:
    BamTagsDict tagsDict(record.tags);
    setTagValue(tagsDict, "BC", barcode);
    writeRecord(bamFileOut, record);
  }
}

#endif /* BARCODE_H_ */


