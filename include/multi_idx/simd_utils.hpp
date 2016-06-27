#pragma once

const __m128i lookup = _mm_setr_epi8(
    /* 0 */ 0, /* 1 */ 1, /* 2 */ 1, /* 3 */ 2,
    /* 4 */ 1, /* 5 */ 2, /* 6 */ 2, /* 7 */ 3,
    /* 8 */ 1, /* 9 */ 2, /* a */ 2, /* b */ 3,
    /* c */ 2, /* d */ 3, /* e */ 3, /* f */ 4
);

const __m128i low_mask = _mm_set1_epi8(0x0f);
const __m128i cleaning_mask = _mm_set1_epi32(0x00ff);

inline __m128i popcount_epi32(const __m128i vec) {
  const __m128i lo      = vec & low_mask;
  const __m128i hi      = _mm_and_si128((_mm_srli_epi16(vec, 4)), low_mask);

  const __m128i popcnt1 = _mm_shuffle_epi8(lookup, lo);
  const __m128i popcnt2 = _mm_shuffle_epi8(lookup, hi);

  __m128i popcounts     = _mm_add_epi8(popcnt1, popcnt2);
  popcounts             = _mm_add_epi8(popcounts, _mm_srli_si128(popcounts, 1));
  popcounts             = _mm_and_si128( _mm_add_epi8(popcounts, _mm_srli_si128(popcounts, 2)), cleaning_mask);
  return popcounts;
}

#define LIKELY(x)   (__builtin_expect((x), 1))
#define UNLIKELY(x) (__builtin_expect((x), 0))