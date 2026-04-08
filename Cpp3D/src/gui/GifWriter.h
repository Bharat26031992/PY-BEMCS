#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <cstring>

namespace BEMCS {

// ============================================================================
// Simple animated GIF writer (no external dependencies)
//
// Usage:
//   GifWriter gif;
//   gif.init(width, height, delay_centiseconds);
//   gif.addFrame(rgbaPixels);   // call for each frame
//   gif.finish("output.gif");
// ============================================================================
class GifWriter {
public:
    void init(int width, int height, int delayCentiseconds = 10) {
        width_ = width;
        height_ = height;
        delay_ = delayCentiseconds;
        frames_.clear();
    }

    // Add a frame from RGBA pixel data (width * height * 4 bytes)
    void addFrame(const uint8_t* rgba) {
        Frame f;
        f.pixels.assign(rgba, rgba + width_ * height_ * 4);
        frames_.push_back(std::move(f));
    }

    // Add a frame from RGB pixel data (width * height * 3 bytes)
    void addFrameRGB(const uint8_t* rgb) {
        Frame f;
        f.pixels.resize(width_ * height_ * 4);
        for (int i = 0; i < width_ * height_; i++) {
            f.pixels[i * 4 + 0] = rgb[i * 3 + 0];
            f.pixels[i * 4 + 1] = rgb[i * 3 + 1];
            f.pixels[i * 4 + 2] = rgb[i * 3 + 2];
            f.pixels[i * 4 + 3] = 255;
        }
        frames_.push_back(std::move(f));
    }

    int frameCount() const { return static_cast<int>(frames_.size()); }

    bool finish(const std::string& filename) {
        if (frames_.empty()) return false;

        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) return false;

        buildGlobalPalette();

        // GIF Header
        out.write("GIF89a", 6);

        // Logical Screen Descriptor
        writeU16(out, width_);
        writeU16(out, height_);
        out.put(static_cast<char>(0xF7)); // GCT flag, 8-bit color, 256 colors
        out.put(0);     // Background color index
        out.put(0);     // Pixel aspect ratio

        // Global Color Table (256 * 3 bytes)
        out.write(reinterpret_cast<const char*>(palette_.data()), 768);

        // Netscape Application Extension (infinite loop)
        out.put(0x21);
        out.put(static_cast<char>(0xFF));
        out.put(11);
        out.write("NETSCAPE2.0", 11);
        out.put(3);
        out.put(1);
        writeU16(out, 0); // 0 = loop forever
        out.put(0);

        // Write each frame
        for (const auto& frame : frames_) {
            // Graphic Control Extension
            out.put(0x21);
            out.put(static_cast<char>(0xF9));
            out.put(4);
            out.put(0x00); // No disposal, no transparency
            writeU16(out, delay_);
            out.put(0);    // No transparent color
            out.put(0);    // Block terminator

            // Image Descriptor
            out.put(0x2C);
            writeU16(out, 0); // Left
            writeU16(out, 0); // Top
            writeU16(out, width_);
            writeU16(out, height_);
            out.put(0); // No local color table

            // Quantize to palette indices
            std::vector<uint8_t> indices(width_ * height_);
            for (int i = 0; i < width_ * height_; i++) {
                indices[i] = findClosestColor(
                    frame.pixels[i * 4 + 0],
                    frame.pixels[i * 4 + 1],
                    frame.pixels[i * 4 + 2]);
            }

            // LZW compress
            writeLZW(out, indices, 8);
        }

        // Trailer
        out.put(0x3B);
        out.close();
        return true;
    }

private:
    struct Frame {
        std::vector<uint8_t> pixels; // RGBA
    };

    int width_ = 0, height_ = 0, delay_ = 10;
    std::vector<Frame> frames_;
    std::vector<uint8_t> palette_; // 768 bytes (256 * RGB)

    void writeU16(std::ofstream& out, int val) {
        out.put(static_cast<char>(val & 0xFF));
        out.put(static_cast<char>((val >> 8) & 0xFF));
    }

    void buildGlobalPalette() {
        palette_.resize(768, 0);
        int idx = 0;
        for (int r = 0; r < 6; r++) {
            for (int g = 0; g < 6; g++) {
                for (int b = 0; b < 6; b++) {
                    if (idx < 256) {
                        palette_[idx * 3 + 0] = static_cast<uint8_t>(r * 51);
                        palette_[idx * 3 + 1] = static_cast<uint8_t>(g * 51);
                        palette_[idx * 3 + 2] = static_cast<uint8_t>(b * 51);
                        idx++;
                    }
                }
            }
        }
        // Fill remaining 40 slots with grayscale ramp
        for (int i = idx; i < 256; i++) {
            uint8_t v = static_cast<uint8_t>((i - idx) * 255 / std::max(1, 255 - idx));
            palette_[i * 3 + 0] = v;
            palette_[i * 3 + 1] = v;
            palette_[i * 3 + 2] = v;
        }
    }

    uint8_t findClosestColor(uint8_t r, uint8_t g, uint8_t b) const {
        int ri = (r + 25) / 51;
        int gi = (g + 25) / 51;
        int bi = (b + 25) / 51;
        ri = std::clamp(ri, 0, 5);
        gi = std::clamp(gi, 0, 5);
        bi = std::clamp(bi, 0, 5);
        return static_cast<uint8_t>(ri * 36 + gi * 6 + bi);
    }

    // Correct LZW compressor using std::unordered_map for the string table
    void writeLZW(std::ofstream& out, const std::vector<uint8_t>& indices,
                  int minCodeSize) {
        out.put(static_cast<char>(minCodeSize));

        const int clearCode = 1 << minCodeSize;
        const int eoiCode = clearCode + 1;

        // Bit-packing state
        uint32_t bitAccum = 0;
        int bitCount = 0;
        std::vector<uint8_t> subBlock;
        subBlock.reserve(260);

        auto writeBits = [&](int code, int numBits) {
            bitAccum |= (static_cast<uint32_t>(code) << bitCount);
            bitCount += numBits;
            while (bitCount >= 8) {
                subBlock.push_back(static_cast<uint8_t>(bitAccum & 0xFF));
                bitAccum >>= 8;
                bitCount -= 8;
                // Flush sub-block when full (max 255 bytes)
                if (subBlock.size() >= 255) {
                    out.put(static_cast<char>(subBlock.size()));
                    out.write(reinterpret_cast<const char*>(subBlock.data()),
                              subBlock.size());
                    subBlock.clear();
                }
            }
        };

        auto flushBits = [&]() {
            if (bitCount > 0) {
                subBlock.push_back(static_cast<uint8_t>(bitAccum & 0xFF));
                bitAccum = 0;
                bitCount = 0;
            }
            if (!subBlock.empty()) {
                out.put(static_cast<char>(subBlock.size()));
                out.write(reinterpret_cast<const char*>(subBlock.data()),
                          subBlock.size());
                subBlock.clear();
            }
            out.put(0); // Block terminator
        };

        // String table: key = (prefix_code << 16) | suffix_byte, value = code
        std::unordered_map<uint32_t, int> table;
        int nextCode;
        int codeSize;

        auto resetTable = [&]() {
            table.clear();
            table.reserve(4096);
            nextCode = eoiCode + 1;
            codeSize = minCodeSize + 1;
        };

        resetTable();
        writeBits(clearCode, codeSize);

        if (indices.empty()) {
            writeBits(eoiCode, codeSize);
            flushBits();
            return;
        }

        int current = indices[0]; // Start with first pixel as initial string

        for (size_t i = 1; i < indices.size(); i++) {
            int pixel = indices[i];
            uint32_t key = (static_cast<uint32_t>(current) << 16) |
                            static_cast<uint32_t>(pixel);

            auto it = table.find(key);
            if (it != table.end()) {
                // String found in table — extend it
                current = it->second;
            } else {
                // String not found — output current code, add new entry
                writeBits(current, codeSize);

                if (nextCode < 4096) {
                    table[key] = nextCode++;
                    // Increase code size when we hit the next power of 2
                    if (nextCode > (1 << codeSize) && codeSize < 12) {
                        codeSize++;
                    }
                } else {
                    // Table full — emit clear code and reset
                    writeBits(clearCode, codeSize);
                    resetTable();
                }

                current = pixel; // Start new string with current pixel
            }
        }

        // Output the last code
        writeBits(current, codeSize);
        writeBits(eoiCode, codeSize);
        flushBits();
    }
};

} // namespace BEMCS
