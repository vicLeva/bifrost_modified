#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "fastq.hpp"

FastqFile::FastqFile(const vector<string> files) : kseq(NULL), fnames(files) {

    fnit = fnames.begin();
    file_no = 0;
    fp = gzopen(fnit->c_str(), "r");
    kseq = kseq_init(fp);
}

FastqFile::~FastqFile() {

  close();
}

void FastqFile::reopen() {

    close();

    fnit = fnames.begin();
    fp = gzopen(fnit->c_str(), "r");
    kseq = kseq_init(fp);
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::read_next(char* read, size_t* read_len, string &seq, size_t* seq_len, unsigned int* file_id, char* qual) {

    int r;

    if ((r = kseq_read(kseq)) >= 0) {

        memcpy(read, kseq->name.s, kseq->name.l + 1); // 0-terminated string
        *read_len = kseq->name.l;
        seq.assign(kseq->seq.s);
        *seq_len = kseq->seq.l;

        if (qual != NULL) memcpy(qual, kseq->qual.s, kseq->qual.l + 1); // 0-terminated string

        if (file_id != NULL) *file_id = file_no / 2;
    }
    else if (r == -1) {

        open_next();

        if (fnit != fnames.end()) return read_next(read, read_len, seq, seq_len, file_id, qual);
        else return -1;
    }

    return r;
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::read_next(string &seq, size_t& id) {

    int r;

    if ((r = kseq_read(kseq)) >= 0) seq.assign(kseq->seq.s);
    else if (r == -1) {

        open_next();

        if (fnit != fnames.end()){

            id = file_no;

            return read_next(seq, id);
        }
    }

    return r;
}

int FastqFile::read_next() {

    int r = kseq_read(kseq);

    if (r == -1) {

        open_next();

        if (fnit != fnames.end()) return read_next();
        return -1;
    }

    return r;
}

vector<string>::const_iterator FastqFile::open_next() {

    if (fnit != fnames.end()) {
        // close current file
        kseq_destroy(kseq);
        gzclose(fp);

        kseq = NULL;

        // get next file
        ++fnit;
        ++file_no;

        if (fnit != fnames.end()) {

            fp = gzopen(fnit->c_str(), "r");
            kseq = kseq_init(fp);
        }
    }

    return fnit;
}

void FastqFile::close() {

    if (kseq != NULL) {

        kseq_destroy(kseq);
        gzclose(fp);

        fnit = fnames.end();
        kseq = NULL;
    }
}
