#include <vector>
#include <iostream>
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
        for (unsigned i = 0; i < (w * h * 3); i++)
        {
            data.push_back((unsigned char)(0));
        }
    }

    PNG(const char *path)
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
        set_pixel(index, r, g, b);
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
        std::vector<unsigned char> png;
        std::cout << "data size: " << data.size() << std::endl;
        std::cout << "width * height * 3: " << width * height * 3 << std::endl;
        unsigned error = lodepng::encode(png, data, width, height, state);
        if (!error)
            error = lodepng::save_file(png, "textures/test_save.png");
        if (error)
            std::cout << "Error saving PNG: " << lodepng_error_text(error) << std::endl;
    };
};