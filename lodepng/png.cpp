#include <vector>
#include "lodepng/lodepng.h"

class PNG
{
public:
    PNG(unsigned int w, unsigned int h)
    {
        width = w;
        height = h;
        data.reserve(sizeof(unsigned char) * w * h * 4);
    }

    PNG(const char* path)
    {
        lodepng::decode(data, width, height, path);
    }
    
    unsigned int width;
    unsigned int height;

    static std::vector<unsigned char> data;

    void set_pixel(unsigned char x,
                   unsigned char y,
                   unsigned char r,
                   unsigned char g,
                   unsigned char b,
                   unsigned char a)
    {
        unsigned long index = (x + (width * y)) * 4;
        data[index] = r;
        data[index + 1] = g;
        data[index + 2] = b;
        data[index + 3] = a;
    };

    void save(const char *path)
    {
        lodepng::encode(path, data, width, height);
    };
};