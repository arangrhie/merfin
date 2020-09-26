TARGET   := merfin
SOURCES  := merfin.C vcf.C varMer.C kmetric.C

SRC_INCDIRS  := . ../utility/src/utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -lmerfin
TGT_PREREQS := libmerfin.a

SUBMAKEFILES :=
