#ifndef gamgee__variant_header__guard
#define gamgee__variant_header__guard

#include "htslib/vcf.h"

#include <memory>
#include <string>
#include <vector>

namespace gamgee {
  
/**
 * @brief Utility class to hold a variant header
 *
 * It can be used to read headers from a VCF/BCF file, but to create one from scratch 
 * you want to use the VariantHeaderBuilder
 */
class VariantHeader {
 public:
  VariantHeader() = default;                                                             ///< @brief initializes a null VariantHeader @warning if you need to create a VariantHeader from scratch, use the builder instead
  explicit VariantHeader(const std::shared_ptr<bcf_hdr_t>& header) : m_header{header} {} ///< @brief creates a VariantHeader given htslib object. @note used by all iterators
  VariantHeader(const VariantHeader& other);                                             ///< @brief makes a deep copy of a VariantHeader. Shared pointers maintain state to all other associated objects correctly.
  VariantHeader(VariantHeader&& other) noexcept;                                         ///< @brief moves VariantHeader accordingly. Shared pointers maintain state to all other associated objects correctly.
  VariantHeader& operator=(const VariantHeader& other);                                  ///< @brief deep copy assignment of a VariantHeader. Shared pointers maintain state to all other associated objects correctly.
  VariantHeader& operator=(VariantHeader&& other) noexcept;                              ///< @brief move assignment of a VariantHeader. Shared pointers maintain state to all other associated objects correctly.
  ~VariantHeader() = default;

  /**
   * @brief returns the number of samples in the header 
   * @note much faster than getting the actual list of samples
   */
  uint32_t n_samples() const { return uint32_t(m_header->n[BCF_DT_SAMPLE]); }  

  std::vector<std::string> filters() const;           ///< @brief builds a vector with the filters
  std::vector<std::string> samples() const;           ///< @brief builds a vector with the names of the samples
  std::vector<std::string> chromosomes() const;       ///< @brief builds a vector with the contigs
  std::vector<std::string> shared_fields() const;     ///< @brief builds a vector with the info fields
  std::vector<std::string> individual_fields() const; ///< @brief builds a vector with the format fields
  bool has_individual_field(const std::string) const; ///< @brief checks if the format field has the given field

  void advanced_merge_header(const VariantHeader& other) { bcf_hdr_combine(other.m_header.get(), m_header.get()); }

  /**
   * @brief looks up the index of a particular filter, shared or individual field tag, enabling subsequent O(1) random-access lookups for that field throughout the iteration. 
   * @return missing_values::int32_t if the tag is not present in the header (you can use missing() on the return value to check)
   * @note prefer this to looking up tag names during the iteration if you are looking for shared fields multiple times. 
   * @note if multiple fields (e.g. shared and individual) have the same tag (e.g. "DP"), they will also have the same index internally, so this function will do the right thing. The accessors for individual and shared field will know how to use the index to retrieve the correct field.
   */
  int32_t field_index(const std::string& tag) const;                                      

 private:
  std::shared_ptr<bcf_hdr_t> m_header;
  
  friend class VariantWriter;
};

}

#endif // gamgee__variant_header__guard
