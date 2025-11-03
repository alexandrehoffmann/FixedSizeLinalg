#ifndef FSLINALG_CRTP_BASE_HPP
#define FSLINALG_CRTP_BASE_HPP

namespace FSLinalg
{
	
template<class Derived>
class CRTPBase
{
public:
	      Derived& derived()       { return static_cast<      Derived&>(*this); }
	const Derived& derived() const { return static_cast<const Derived&>(*this); }
};
	
} // namespace FSLinalg

#endif // FSLINALG_CRTP_BASE_HPP
