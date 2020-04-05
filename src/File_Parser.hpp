#ifndef BIFROST_FILE_PARSER_HPP
#define BIFROST_FILE_PARSER_HPP

#include <sstream>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "FASTX_Parser.hpp"
#include "GFA_Parser.hpp"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif

class FileParser {

    public:

        FileParser(const vector<string>& filenames) :   files_it(0), files_fastx_it(0), files_gfa_it(0),
                                                        reading_fastx(false), invalid(false) {

            if (filenames.size() == 0) {

                cerr << "FileParser::FileParser(): Missing input files" << endl;
                invalid = true;
            }
            else {

                struct stat stFileInfo;

                files = filenames;

                for (const auto& s : files) {

                    const int intStat = stat(s.c_str(), &stFileInfo);

                    if (intStat != 0) {

                        cerr << "FileParser::FileParser(): File not found: " << s << endl;
                        invalid = true;
                    }
                    else {

                        const int format = FileParser::getFileFormat(s.c_str());

                        if (format == -1){

                            cerr << "FileParser::FileParser(): Input file " << s << " does not exist, is ill-formed or is not in FASTA/FASTQ/GFA format." << endl;

                            invalid = true;
                        }
                        else if (format == 0) files_fastx.push_back(s); // FASTA
                        else if (format == 1) files_fastx.push_back(s); // FASTQ
                        else if (format == 2) files_gfa.push_back(s); // GFA
                    }
                }
            }

            if (!invalid){

                if (!files_fastx.empty()){

                    ff = FastqFile(files_fastx);
                    reading_fastx = (files[0] == files_fastx[0]);
                }

                if (!files_gfa.empty()){

                    gfap = GFA_Parser(files_gfa);
                    invalid = !gfap.open_read();
                    reading_fastx = (files[0] != files_gfa[0]);
                }
            }
        }

        bool read(string& seq, size_t& file_id){

            if (!invalid){

                bool new_file;

                if (reading_fastx){

                    const int ret = ff.read_next(seq, files_fastx_it, new_file);

                    if (new_file || (ret == -1)){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = ((ret != -1) && !files_fastx.empty() && (files[files_it] == files_fastx[files_fastx_it]));

                            return read(seq, file_id); // We read the next line of the file
                        }
                    }
                }
                else {
                    // Read first line of next GFA file, skip edge lines
                    GFA_Parser::GFA_line r = gfap.read(files_gfa_it, new_file, true);

                    if (new_file || ((r.first == nullptr) && (r.second == nullptr))){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = !files_fastx.empty() && (files[files_it] == files_fastx[files_fastx_it]);

                            return read(seq, file_id); // We read the next line of the file
                        }
                    }

                    if (r.first != nullptr) seq = r.first->seq;
                }

                file_id = files_it;
            }

            return !invalid;
        }

        bool read(stringstream& ss, size_t& file_id){

            if (!invalid){

                bool new_file;

                if (reading_fastx){

                    const int ret = ff.read_next(ss, files_fastx_it, new_file);

                    if (new_file || (ret == -1)){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = ((ret != -1) && !files_fastx.empty() && (files[files_it] == files_fastx[files_fastx_it]));

                            return read(ss, file_id); // We read the next line of the file
                        }
                    }
                }
                else {
                    // Read first line of next GFA file, skip edge lines
                    GFA_Parser::GFA_line r = gfap.read(files_gfa_it, new_file, true);

                    if (new_file || ((r.first == nullptr) && (r.second == nullptr))){ // Need to open next file

                        invalid = ((files_it + 1) >= files.size()); // Invalidate this object if no more file to read

                        if (!invalid){ // If still some files to read

                            ++files_it; //Increment iterator to next file to read

                            // Check if next file to read is FASTA/FASTQ format
                            reading_fastx = !files_fastx.empty() && (files[files_it] == files_fastx[files_fastx_it]);

                            return read(ss, file_id); // We read the next line of the file
                        }
                    }

                    if (r.first != nullptr) ss << r.first->seq;
                }

                file_id = files_it;
            }

            return !invalid;
        }

        const char* getNameString() const {

            if (invalid || !reading_fastx) return nullptr;
            return ff.get_kseq()->name.s;
        }

        const char* getQualityScoreString() const {

            if (invalid || !reading_fastx) return nullptr;
            return ff.get_kseq()->qual.s;
        }

        void close(){

            ff.close();
            gfap.close();
        }

        // Returns 0 for FASTA, 1 for FASTQ, 2 for GFA, -1 otherwise
        static int getFileFormat(const char* filename) {

            char gfa_header[] = "H\tVN:Z:";

            const size_t sz_gfa_header = strlen(static_cast<char*>(gfa_header));
            const size_t sz_buffer = 16384;

            char buffer[sz_buffer];

            int ret_v = -1;

            gzFile fp = gzopen(filename, "r");

            if (fp == Z_NULL) return ret_v; // Couldn't open file, return undetermined file format

            const int l_read = gzread(fp, buffer, sz_buffer - 1);

            buffer[l_read] = '\0';

            if (l_read != 0) { // Must contain at least one character, otherwise return undertermined file format

                if ((buffer[0] == '>') || (buffer[0] == '@')) { // File is FASTA or FASTQ format

                    gzrewind(fp); // Put back file cursor to beginning;

                    kseq_t* kseq = kseq_init(fp); // Initialize kseq

                    if (kseq != NULL){

                        const int r = kseq_read(kseq); // Read first record of FASTA or FASTQ file

                        // If file contains at least 1 record with a name and sequence, returns 0: FASTA, 1: FASTQ
                        if ((r >= 0) && (kseq->name.l >= 1) && (kseq->seq.l >= 1)) ret_v = static_cast<int>(kseq->qual.l == kseq->seq.l);

                        kseq_destroy(kseq);
                    }
                }
                else if (strncmp(buffer, static_cast<char*>(gfa_header), sz_gfa_header) == 0) ret_v = 2; // GFA
            }

            gzclose(fp);

            return ret_v;
        }

    private:

        bool invalid;
        bool reading_fastx;

        size_t files_it;
        size_t files_fastx_it;
        size_t files_gfa_it;

        vector<string> files;
        vector<string> files_fastx;
        vector<string> files_gfa;

        FastqFile ff;
        GFA_Parser gfap;
};

#endif
