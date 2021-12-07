#ifndef FILE_RANDOM_HPP
#define FILE_RANDOM_HPP
#include <cstdint>
#include <limits>

struct U32Random {
   std::uint64_t state{0};
   std::uint64_t inc{seed()};

   using result_type = std::uint32_t;

   [[nodiscard]]
   constexpr result_type operator()() noexcept {
      return next();
   }

   [[nodiscard]]
   static constexpr result_type min() noexcept {
      return std::numeric_limits<result_type>::min();
   }
   [[nodiscard]]
   static constexpr result_type max() noexcept {
      return std::numeric_limits<result_type>::max();
   }
private:
   [[nodiscard]]
   constexpr result_type next() noexcept {
      const auto old_state = state;
      const auto xor_shifted = ((old_state >> 18u) ^ old_state) >> 27u;
      const auto rot = old_state >> 59u;
      state = old_state * 6364136223846793005ULL + (inc | 1);
      return static_cast<std::uint32_t>((xor_shifted >> rot) | (xor_shifted << ((-rot) & 31)));
   }
   [[nodiscard]]
   static constexpr std::uint64_t seed() noexcept {
      std::uint64_t shifted = 0;
      for (const auto c : __TIME__) {
         shifted <<= 8;
         shifted  |= c;
      }
      return shifted;
   }
};

struct F32Random {
   U32Random rand{};

   using result_type = float;

   [[nodiscard]]
   constexpr result_type operator()() noexcept {
      return ((result_type)rand() / (result_type)U32Random::max());
   }

   [[nodiscard]]
   static constexpr result_type min() noexcept {
      return 0.0f;
   }

   [[nodiscard]]
   static constexpr result_type max() noexcept {
      return 1.0f;
   }
};


#endif /* FILE_RANDOM_HPP */
