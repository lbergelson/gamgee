#include "reference_block_splitting_variant_iterator.h"

namespace gamgee {

ReferenceBlockSplittingVariantIterator::ReferenceBlockSplittingVariantIterator() :
  MultipleVariantIterator {},
  m_pending_variants {},
  m_split_variants {}
{}

ReferenceBlockSplittingVariantIterator::ReferenceBlockSplittingVariantIterator(const std::vector<vcfFile*> variant_files, const std::shared_ptr<bcf_hdr_t> variant_header) :
  MultipleVariantIterator {variant_files, variant_header},
  m_pending_variants {},
  m_split_variants {}
{
  m_split_variants.reserve(variant_files.size());
  fetch_next_split_vector();
}

ReferenceBlockSplittingVariantIterator::ReferenceBlockSplittingVariantIterator(ReferenceBlockSplittingVariantIterator&& original) :
  MultipleVariantIterator {std::move(original)},
  m_pending_variants {std::move(original.m_pending_variants)},
  m_split_variants {std::move(original.m_split_variants)}
{}

std::vector<Variant>& ReferenceBlockSplittingVariantIterator::operator*() {
  return m_split_variants;
}

std::vector<Variant>& ReferenceBlockSplittingVariantIterator::operator++() {
  fetch_next_split_vector();
  return m_split_variants;
}

// NOTE: this method does the minimal work necessary to determine that we have reached the end of iteration
// it is NOT a valid general-purpose inequality method
bool ReferenceBlockSplittingVariantIterator::operator!=(const ReferenceBlockSplittingVariantIterator& rhs) {
  return !(m_pending_variants.empty() && rhs.m_pending_variants.empty()
      && m_split_variants.empty() && rhs.m_split_variants.empty());
}

void ReferenceBlockSplittingVariantIterator::populate_pending () {
  for (const auto& variant : MultipleVariantIterator::operator*()) {
    m_pending_min_end = std::min(m_pending_min_end, variant.alignment_stop());
    m_pending_variants.push_front(std::move(variant));
  }

  MultipleVariantIterator::operator++();
}

void ReferenceBlockSplittingVariantIterator::populate_split_variants () {
  // m_pending_chrom can only change when pending is empty
  auto new_pending_start = -1;
  auto new_pending_end = UINT_MAX;

  // advance iter by erase() or ++ depending on if
  for (auto iter = m_pending_variants.begin(); iter != m_pending_variants.end(); ) {
    auto& variant = *iter;
    auto var_end = variant.alignment_stop();
    if (var_end == m_pending_min_end) {
      m_split_variants.push_back(std::move(variant));
      iter = m_pending_variants.erase(iter);
    }
    else {
      auto split_variant = variant;
      split_variant.set_alignment_stop(m_pending_min_end);
      m_split_variants.push_back(std::move(split_variant));

      new_pending_start = m_pending_min_end + 1;
      new_pending_end = std::min(new_pending_end, var_end);
      variant.set_alignment_start(new_pending_start);
      variant.set_alignment_stop(var_end);      // stop is internally an offset to start, so we need to reset it after updating stop
      iter++;
    }
  }

  if (new_pending_start != -1) {
    m_pending_start = new_pending_start;
    m_pending_min_end = new_pending_end;
  }

}

void ReferenceBlockSplittingVariantIterator::fetch_next_split_vector() {
  m_split_variants.clear();
  // run until we have split variants or the incoming pre-split vector is empty, indicating iteration is done
  // if the incoming vector is empty, we also end with an empty vector to indicate our iteration is done
  while (m_split_variants.empty() && ! MultipleVariantIterator::operator*().empty()) {
    auto& incoming = MultipleVariantIterator::operator*();
    // when the pending list is empty, fill from the incoming vector and consume it
    if (m_pending_variants.empty()) {
      m_pending_chrom = incoming[0].chromosome();
      m_pending_start = incoming[0].alignment_start();
      m_pending_min_end = UINT_MAX;
      populate_pending();
    }
    else {
      // if the incoming vector has the same start location as pending, add it to pending and consume it
      if (!incoming.empty()
          && incoming[0].chromosome() == m_pending_chrom
          && incoming[0].alignment_start() == m_pending_start)
        populate_pending();
    }

    // refresh incoming because it may have advanced
    incoming = MultipleVariantIterator::operator*();

    // if the pending Variants' end is after the incoming vector's start, then their end is too high
    if (!incoming.empty() && incoming[0].chromosome() == m_pending_chrom)
      m_pending_min_end = std::min(m_pending_min_end, incoming[0].alignment_start() - 1);

    populate_split_variants();
  }
}

}

