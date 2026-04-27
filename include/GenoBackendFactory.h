/*
   GCTA: factory functions that create concrete GenoBackend instances.

   Each factory lives in its own translation unit (BedBackend.cpp, etc.)
   so that the concrete class definitions never need to be visible to
   Geno.cpp — only the factory signatures do.

   Geno& must already have its format-independent state populated
   (sampleKeepIndex, rawSampleCT, geno_files, marker, …) before a factory
   is called.  The factory constructor performs the format-specific
   pre-logic (allocating masks, computing strides, etc.) and its destructor
   performs the matching cleanup.
*/
#pragma once

#include <memory>
#include "GenoBackend.h"

// Forward declaration only — concrete class definitions are in the .cpp files.
class Geno;

/// Create a BED-format backend.  Also handles the PGEN-aliased-to-BED case.
std::unique_ptr<GenoBackend> makeBedBackend(Geno &geno);

/// Create a BGEN-format backend.
std::unique_ptr<GenoBackend> makeBgenBackend(Geno &geno);

/// Create a PGEN-format backend (proper PGEN via ReadDosage).
std::unique_ptr<GenoBackend> makePgenBackend(Geno &geno);

/// Create a VCF/BCF-format backend (htslib-based, random access via CSI/TBI index).
std::unique_ptr<GenoBackend> makeVcfBackend(Geno &geno);
