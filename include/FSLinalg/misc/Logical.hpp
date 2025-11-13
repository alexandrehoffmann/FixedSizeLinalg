#ifndef FSLINALG_MISC_LOGICAL_HPP
#define FSLINALG_MISC_LOGICAL_HPP

namespace FSLinalg
{
	constexpr bool implies(const bool p, const bool q) { return (not p) or q; }
}

#endif // FSLINALG_MISC_LOGICAL_HPP
