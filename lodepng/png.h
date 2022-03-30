#include <vector>
#include "lodepng.h"

class PNG
{
public:
    unsigned int width;
    unsigned int height;

    std::vector<unsigned char> data;

    PNG(unsigned int w, unsigned int h)
    {
        width = w;
        height = h;
        data.reserve(sizeof(unsigned char) * w * h * 3);
    }

    PNG(const char* path)
    {
        lodepng::decode(data, width, height, path);
    }

    void set_pixel(unsigned int x,
                   unsigned int y,
                   unsigned char r,
                   unsigned char g,
                   unsigned char b)
    {
        unsigned long index = (x + (width * y)) * 3;
        data[index] = r;
        data[index + 1] = g;
        data[index + 2] = b;
    }

    void set_pixel(unsigned long i,
                   unsigned char r,
                   unsigned char g,
                   unsigned char b)
    {
        data[i] = r;
        data[i + 1] = g;
        data[i + 2] = b;
    };

    void save(const char *path)
    {
        lodepng::State state;
        state.info_raw.colortype = LCT_RGB;
        lodepng::encode(path, data, width, height);
    };
};