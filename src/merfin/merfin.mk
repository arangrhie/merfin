TARGET   := merfin
SOURCES  := merfin.C \
            merfin-globals.C \
            merfin-histogram.C \
            merfin-dump.C \
            merfin-variants.C \
            merfin-completeness.C \
            vcf.C \
            varMer.C

SRC_INCDIRS  := . ../utility/src/utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lmerfin
TGT_PREREQS := libmerfin.a

SUBMAKEFILES :=
