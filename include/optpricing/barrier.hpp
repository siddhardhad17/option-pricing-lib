#pragma once

#include "optpricing/types.hpp"

namespace optpricing {

// Closed-form price of a continuously monitored barrier option with zero rebate
// (Reiner & Rubinstein, 1991; see Haug, "The Complete Guide to Option Pricing Formulas").
//
// If the barrier has already been breached at inception, the knock-out value is 0 and the
// knock-in value equals the vanilla price.
double barrier_closed_form(const Market& m, const BarrierOption& o, double vol);

}  // namespace optpricing
