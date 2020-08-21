#pragma once

class MathUtil {
    public:
        static float interp(float x, float x0, float x1, float y0, float y1)
        {
            return lerp(y0, y1, alpha(x0, x1, x));
        }

        static float alpha(float x0, float x1, float x)
        {
            return (x - x0) / (x1 - x0);
        }

        static float lerp(float a, float b, float t)
        {
            return a + t * (b - a);
        }
};