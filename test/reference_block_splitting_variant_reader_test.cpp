#include "multiple_variant_reader.h"
#include "reference_block_splitting_variant_iterator.h"

#include <boost/test/unit_test.hpp>

using namespace std;
using namespace gamgee;

BOOST_AUTO_TEST_CASE( split_reference_blocks )
{
  const auto test_files = vector<string>{"testdata/ref_block/test1.vcf", "testdata/ref_block/test2.vcf", "testdata/ref_block/test3.vcf",
    "testdata/ref_block/test4.vcf", "testdata/ref_block/test5.vcf"};
  const auto truth_contigs = vector<uint32_t>{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  const auto truth_block_starts = vector<uint32_t>{1, 20, 41, 42, 50, 61, 81, 101, 102, 103, 151, 201, 300, 401, 600, 650, 701, 10001000};
  const auto truth_block_stops = vector<uint32_t>{19, 40, 41, 49, 60, 80, 100, 101, 102, 150, 200, 299, 400, 500, 649, 700, 750, 10001001};
  // note: final block is not a ref block, but the block ends at start + 1 because the reference allele length is 2
  GVCFReader reader{test_files, false};
  auto position_counter = 0u;
  for (const auto& vec : reader) {
    for (const auto& record : vec) {
      BOOST_CHECK_EQUAL(record.chromosome(), truth_contigs[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_start(), truth_block_starts[position_counter]);
      BOOST_CHECK_EQUAL(record.alignment_stop(), truth_block_stops[position_counter]);
    }
    ++position_counter;
  }
  BOOST_CHECK_EQUAL(position_counter, 18u);
}
