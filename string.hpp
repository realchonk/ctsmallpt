#ifndef FILE_STRING_HPP
#define FILE_STRING_HPP
#include <string_view>

// constexpr string

constexpr auto copy(auto* dest, const auto* src, const std::size_t len) {
   for (std::size_t i = 0; i < len; ++i) {
      dest[i] = src[i];
   }
   return dest;
}

class String {
private:
   std::size_t len{}, cap{};
   char* buffer{};
public:
   constexpr String() noexcept = default;

   constexpr String(const String& s)
      : len(s.len), cap(s.cap), buffer(new char[cap]) {
      copy(buffer, s.buffer, s.len + 1);
   }

   constexpr String(String&& s) noexcept
      : len(s.len), cap(s.cap), buffer(s.buffer) {
      s = String();
   }

   constexpr String(std::string_view s)
      : len(s.size()), cap(s.size() + 1), buffer(new char[cap]) {
      copy(buffer, s.data(), len);
      buffer[len] = '\0';
   }

   constexpr ~String() noexcept {
      delete[] buffer;
      len = cap = 0;
      buffer = nullptr;
   }

   // Assignment operators

   constexpr String& operator=(const String& s) {
      reserve(s.len + 1, false);
      len = s.len;
      copy(buffer, s.buffer, len + 1);
      return *this;
   }
   constexpr String& operator=(String&& s) noexcept {
      delete[] buffer;
      len = s.len;
      cap = s.cap;
      buffer = s.buffer;
      s.buffer = nullptr;
      s.len = s.cap = 0;
      return *this;
   }
   constexpr String& operator=(std::string_view s) {
      len = s.size();
      reserve(len + 1, false);
      copy(buffer, s.data(), s.size());
      buffer[len] = '\0';
      return *this;
   }

   // Append operators

   constexpr String& operator+=(const String& s) {
      reserve(len + s.len + 1, true);
      copy(buffer + len, s.buffer, s.len);
      len += s.len;
      buffer[len] = '\0';
      return *this;
   }

   constexpr String& operator+=(std::string_view s) {
      reserve(len + s.size() + 1, true);
      copy(buffer + len, s.data(), s.size());
      len += s.size();
      buffer[len] = '\0';
      return *this;
   }

   constexpr String& operator+=(char ch) {
      reserve(len + 2, true);
      buffer[len    ] = ch;
      buffer[len + 1] = '\0';
      len += 1;
      return *this;
   }

   // Index operator

   constexpr char& operator[](std::size_t i) {
      if (i < len) {
         return buffer[i];
      } else {
         throw std::runtime_error("Out of bounds");
      }
   }
   constexpr char operator[](std::size_t i) const {
      if (i < len) {
         return buffer[i];
      } else {
         throw std::runtime_error("Out of bounds");
      }
   }

   [[nodiscard]]
   constexpr operator std::string_view() const noexcept {
      return { buffer, buffer + len };
   }

   [[nodiscard]]
   constexpr auto size() const noexcept { return len; }

   [[nodiscard]]
   constexpr bool empty() const noexcept { return !len; }

   [[nodiscard]]
   constexpr auto capacity() const noexcept { return cap; }

   [[nodiscard]]
   constexpr auto data() const noexcept { return buffer; }

   [[nodiscard]]
   constexpr auto c_str() const noexcept { return buffer; }

   [[nodiscard]]
   constexpr auto begin() noexcept { return buffer; }

   [[nodiscard]]
   constexpr auto begin() const noexcept { return buffer; }

   [[nodiscard]]
   constexpr auto cbegin() const noexcept { return buffer; }

   [[nodiscard]]
   constexpr auto end() noexcept { return buffer + len; }

   [[nodiscard]]
   constexpr auto end() const noexcept { return buffer + len; }

   [[nodiscard]]
   constexpr auto cend() const noexcept { return buffer + len; }


   constexpr void clear() noexcept {
      buffer[0] = '\0';
      len = 0;
   }

   constexpr void reserve(std::size_t num, bool do_copy = true) {
      if (cap < num) {
         cap = num;
         auto old_buffer = buffer;
         buffer = new char[cap];
         if (do_copy)
            copy(buffer, old_buffer, len);
         delete[] old_buffer;
      }
   }

   constexpr String& reverse() noexcept {
      if (empty())
         return *this;
      auto a = begin();
      auto b = end();
      while (a < b) {
         --b;
         std::swap(*a, *b);
         ++a;
      }
      return *this;
   }

   static constexpr String parse_int(std::integral auto val) {
      String str{};
      auto negative = false;
      if (val < 0) {
         negative = true;
         val = -val;
      }
      while (val != 0) {
         str += (char)(val % 10 + '0');
         val /= 10;
      }
      if (negative)
         str += '-';
      str.reverse();
      return str;
   }

};

#endif /* FILE_STRING_HPP */
