#include "sam.h"
#include "cigar.h"
#include "utils/hts_memory.h"
#include "htslib/sam.h"
#include "missing.h"

#include <iostream>
#include <string>

using namespace std;

namespace gamgee {

constexpr auto MATE_CIGAR_TAG = "MC";

Sam::Sam(const std::shared_ptr<bam_hdr_t>& header, const std::shared_ptr<bam1_t>& body) noexcept :
  m_header {header},
  m_body {body}
{}

Sam::Sam(const Sam& other) :
  m_header { other.m_header },
  m_body { utils::make_shared_sam(utils::sam_deep_copy(other.m_body.get())) }
{}

Sam& Sam::operator=(const Sam& other) {
  if ( &other == this )  
    return *this;
  m_header = other.m_header;      ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  m_body = utils::make_shared_sam(utils::sam_deep_copy(other.m_body.get()));     ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  return *this;
}

uint32_t Sam::mate_alignment_stop(const SamTag<string>& mate_cigar_tag) const {
  auto result = mate_alignment_start();
  stringstream cigar_stream {mate_cigar_tag.value()};
  auto has_reference_bases = false;
  while (cigar_stream.peek() != std::char_traits<char>::eof()) {
    const auto element = Cigar::parse_next_cigar_element(cigar_stream);
    if (Cigar::consumes_reference_bases(Cigar::cigar_op(element))) {
      result +=  Cigar::cigar_oplen(element);
      has_reference_bases = true;
    }
  }
  return has_reference_bases ? result-1 : result; // we want the last base to be 1-based and **inclusive** but we only need to deduct 1 if we had any reference consuming operators in the read.
}

uint32_t Sam::mate_alignment_stop() const {
  const auto mate_cigar = string_tag(MATE_CIGAR_TAG);
  if (missing(mate_cigar))
    throw std::invalid_argument{string{"Cannot find the mate alignment stop on a record without the tag: "} + MATE_CIGAR_TAG};
  return mate_alignment_stop(mate_cigar);
}

uint32_t Sam::unclipped_start() const {
  auto pos = alignment_start();
  const auto* cigar = bam_get_cigar(m_body.get());
  for (auto i = 0u; i != m_body->core.n_cigar; ++i) {
    const auto op = bam_cigar_op(cigar[i]);
    if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
      pos -= bam_cigar_oplen(cigar[i]);
    else
      break;
  }
  return pos;
}

uint32_t Sam::unclipped_stop() const {
  auto pos = alignment_stop();
  const auto* cigar = bam_get_cigar(m_body.get());
  for (auto i = int{m_body->core.n_cigar - 1}; i >= 0; --i) {
    const auto op = bam_cigar_op(cigar[i]);
    if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)
      pos += bam_cigar_oplen(cigar[i]);
    else
      break;
  }
  return pos;
}

uint32_t Sam::mate_unclipped_start() const {
  const auto mate_cigar_tag = string_tag(MATE_CIGAR_TAG);
  if (missing(mate_cigar_tag))
    throw std::invalid_argument{string{"Cannot find the mate unclipped start on a record without the tag: "} + MATE_CIGAR_TAG};
  return mate_unclipped_start(mate_cigar_tag);
}

uint32_t Sam::mate_unclipped_start(const SamTag<string>& mate_cigar_tag) const {
  auto result = mate_alignment_start();
  stringstream cigar_stream {mate_cigar_tag.value()};
  while (cigar_stream.peek() != std::char_traits<char>::eof()) {
    const auto element = Cigar::parse_next_cigar_element(cigar_stream);
    const auto op = Cigar::cigar_op(element);
    if (op != CigarOperator::S && op != CigarOperator::H)
      break;
    result -= Cigar::cigar_oplen(element);
  }
  return result;
}

uint32_t Sam::mate_unclipped_stop() const {
  const auto mate_cigar_tag = string_tag(MATE_CIGAR_TAG);
  if (missing(mate_cigar_tag))
    throw std::invalid_argument{string{"Cannot find the mate unclipped stop on a record without the tag: "} + MATE_CIGAR_TAG};
  return mate_unclipped_stop(mate_cigar_tag);
}

uint32_t Sam::mate_unclipped_stop(const SamTag<string>& mate_cigar_tag) const {
  auto result = mate_alignment_start();
  stringstream cigar_stream {mate_cigar_tag.value()};
  auto end_tail = false;
  auto has_reference_bases = false;
  while (cigar_stream.peek() != std::char_traits<char>::eof()) {
    const auto element = Cigar::parse_next_cigar_element(cigar_stream);
    const auto op = Cigar::cigar_op(element);
    if (!end_tail && (op == CigarOperator::S || op == CigarOperator::H))
      continue;
    if (!end_tail)
      end_tail = true;
    if (Cigar::consumes_reference_bases(op) || op == CigarOperator::S || op == CigarOperator::H) {
      result += Cigar::cigar_oplen(element);
      has_reference_bases = true;
    }
  }
  return has_reference_bases ? result - 1 : result; // we only want to subtract 1 in case we actually had bases in the read that modified the alignment start. Because we are adding the length of the operators, we will end up pointing to one past the unclipped stop of the read. We want to be **inclusive**. 
}



/**
 * @brief retrieve an integer-valued tag by name
 *
 * @note returns a SamTag with missing() == true if the read has no tag by this name
 */
SamTag<int32_t> Sam::integer_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  return aux_ptr == nullptr ? SamTag<int32_t>(tag_name, 0, true) :
                              SamTag<int32_t>(tag_name, bam_aux2i(aux_ptr));

  // TODO: htslib bam_aux2i() returns 0 if the tag is not of integer type,
  //       but 0 is a valid value for a correctly-typed tag. This should be patched.
}

/**
 * @brief retrieve a double/float-valued tag by name
 *
 * @note returns a SamTag with missing() == true if the read has no tag by this name
 */
SamTag<double> Sam::double_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  return aux_ptr == nullptr ? SamTag<double>(tag_name, 0.0, true) :
                              SamTag<double>(tag_name, bam_aux2f(aux_ptr));

  // TODO: htslib bam_aux2f() returns 0.0 if the tag is not of double/float type,
  //       but 0.0 is a valid value for a correctly-typed tag. This should be patched.
}

/**
 * @brief retrieve a char-valued tag by name
 *
 * @note returns a SamTag with missing() == true if the read has no tag by this name
 */
SamTag<char> Sam::char_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  if ( aux_ptr == nullptr )  // tag doesn't exist
    return SamTag<char>(tag_name, '\0', true);

  const auto char_val = bam_aux2A(aux_ptr);
  if ( char_val == '\0' )    // tag not of char type
    return SamTag<char>(tag_name, '\0', true);

  return SamTag<char>(tag_name, char_val);
}

/**
 * @brief retrieve a string-valued tag by name
 *
 * @note returns a SamTag with missing() == true if the read has no tag by this name
 */
SamTag<std::string> Sam::string_tag(const std::string& tag_name) const {
  const auto aux_ptr = bam_aux_get(m_body.get(), tag_name.c_str());
  if ( aux_ptr == nullptr )  // tag doesn't exist
    return SamTag<string>(tag_name, "", true);

  const auto str_ptr = bam_aux2Z(aux_ptr);
  if ( str_ptr == nullptr )  // tag not of string type
    return SamTag<string>(tag_name, "", true);

  return SamTag<string>(tag_name, string{str_ptr});
}

}

