package net.minecraft.src;

import java.util.Random;

public class NoiseGeneratorPerlin extends NoiseGenerator {
    private int[] permutations;
    public double xCoord;
    public double yCoord;
    public double zCoord;
    private Noise noise;

    static class NoiseParams {
        float offset;
        float scale;
        Vector3 spread;
        int seed;
        int octaves;
        float persist;
        float lacunarity;
        int flags;

        NoiseParams() {
            offset = 0.0f;
            scale = 1.0f;
            spread = new Vector3(600.0f, 600.0f, 150.0f); // Adjusted for Minetest-like terrain
            seed = 0;
            octaves = 3; // Increased for fractal detail
            persist = 0.6f; // Adjusted for smoother transitions
            lacunarity = 2.0f;
            flags = 2; // Enable NOISE_FLAG_EASED
        }
    }

    static class Vector3 {
        float X, Y, Z;

        Vector3(float x, float y, float z) {
            X = x;
            Y = y;
            Z = z;
        }
    }

    static class Noise {
        private NoiseParams np;
        private int seed;
        private int sx, sy, sz;
        private float[] value_buf;
        private float[] persist_buf;
        private float[] noise_buf;
        private float[] result;

        Noise(NoiseParams np, int seed, int sx, int sy, int sz) {
            this.np = np;
            this.seed = seed;
            this.sx = sx;
            this.sy = sy;
            this.sz = sz;
            allocBuffers();
        }

        void allocBuffers() {
            if (sx < 1) sx = 1;
            if (sy < 1) sy = 1;
            if (sz < 1) sz = 1;

            int bufsize = sx * sy * sz;
            value_buf = new float[bufsize];
            result = new float[bufsize];
            persist_buf = null;
            resizeNoiseBuf(sz > 1);
        }

        void resizeNoiseBuf(boolean is3d) {
            float ofactor = np.lacunarity > 1.0f ? (float) Math.pow(np.lacunarity, np.octaves - 1) : np.lacunarity;
            float num_noise_points_x = sx * ofactor / np.spread.X;
            float num_noise_points_y = sy * ofactor / np.spread.Y;
            float num_noise_points_z = sz * ofactor / np.spread.Z;

            if (num_noise_points_x > 1000000000.0f || num_noise_points_y > 1000000000.0f || num_noise_points_z > 1000000000.0f) {
                throw new RuntimeException("Invalid noise parameters");
            }

            if (np.spread.X / ofactor < 1.0f || np.spread.Y / ofactor < 1.0f || np.spread.Z / ofactor < 1.0f) {
                throw new RuntimeException("Too many octaves");
            }

            int nlx = (int) Math.ceil(num_noise_points_x) + 3;
            int nly = (int) Math.ceil(num_noise_points_y) + 3;
            int nlz = is3d ? (int) Math.ceil(num_noise_points_z) + 3 : 1;

            noise_buf = new float[nlx * nly * nlz];
        }

        float noise2d(int x, int y, int seed, int[] permutations) {
            int n = (1619 * x + 31337 * y + 1013 * seed) & 0x7fffffff;
            n = permutations[n & 255];
            n = (n >> 13) ^ n;
            n = (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff;
            return 1.0f - (float) n / 0x40000000;
        }

        float noise3d(int x, int y, int z, int seed, int[] permutations) {
            int n = (1619 * x + 31337 * y + 52591 * z + 1013 * seed) & 0x7fffffff;
            n = permutations[n & 255];
            n = (n >> 13) ^ n;
            n = (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff;
            return 1.0f - (float) n / 0x40000000;
        }

        float noise2d_value(float x, float y, int seed, boolean eased, int[] permutations) {
            int x0 = (int) Math.floor(x);
            int y0 = (int) Math.floor(y);
            float xl = x - x0;
            float yl = y - y0;

            float v00 = noise2d(x0, y0, seed, permutations);
            float v10 = noise2d(x0 + 1, y0, seed, permutations);
            float v01 = noise2d(x0, y0 + 1, seed, permutations);
            float v11 = noise2d(x0 + 1, y0 + 1, seed, permutations);

            return biLinearInterpolation(v00, v10, v01, v11, xl, yl, eased);
        }

        float noise3d_value(float x, float y, float z, int seed, boolean eased, int[] permutations) {
            int x0 = (int) Math.floor(x);
            int y0 = (int) Math.floor(y);
            int z0 = (int) Math.floor(z);
            float xl = x - x0;
            float yl = y - y0;
            float zl = z - z0;

            float v000 = noise3d(x0, y0, z0, seed, permutations);
            float v100 = noise3d(x0 + 1, y0, z0, seed, permutations);
            float v010 = noise3d(x0, y0 + 1, z0, seed, permutations);
            float v110 = noise3d(x0 + 1, y0 + 1, z0, seed, permutations);
            float v001 = noise3d(x0, y0, z0 + 1, seed, permutations);
            float v101 = noise3d(x0 + 1, y0, z0 + 1, seed, permutations);
            float v011 = noise3d(x0, y0 + 1, z0 + 1, seed, permutations);
            float v111 = noise3d(x0 + 1, y0 + 1, z0 + 1, seed, permutations);

            return triLinearInterpolation(v000, v100, v010, v110, v001, v101, v011, v111, xl, yl, zl, eased);
        }

        float biLinearInterpolation(float v00, float v10, float v01, float v11, float x, float y, boolean eased) {
            if (eased) {
                x = easeCurve(x);
                y = easeCurve(y);
            }
            float u = linearInterpolation(v00, v10, x);
            float v = linearInterpolation(v01, v11, x);
            return linearInterpolation(u, v, y);
        }

        float triLinearInterpolation(float v000, float v100, float v010, float v110,
                                     float v001, float v101, float v011, float v111,
                                     float x, float y, float z, boolean eased) {
            if (eased) {
                x = easeCurve(x);
                y = easeCurve(y);
                z = easeCurve(z);
            }
            float u = biLinearInterpolation(v000, v100, v010, v110, x, y, eased);
            float v = biLinearInterpolation(v001, v101, v011, v111, x, y, eased);
            return linearInterpolation(u, v, z);
        }

        float linearInterpolation(float v0, float v1, float t) {
            return v0 + (v1 - v0) * t;
        }

        float easeCurve(float t) {
            return t * t * (3.0f - 2.0f * t);
        }

        float[] noiseMap2D(float x, float y, float[] persistence_map, int[] permutations) {
            float f = 1.0f, g = 1.0f;
            int bufsize = sx * sy;

            x /= np.spread.X;
            y /= np.spread.Y;

            for (int i = 0; i < bufsize; i++) {
                result[i] = 0.0f;
            }

            if (persistence_map != null) {
                if (persist_buf == null) {
                    persist_buf = new float[bufsize];
                }
                for (int i = 0; i < bufsize; i++) {
                    persist_buf[i] = 1.0f;
                }
            }

            for (int oct = 0; oct < np.octaves; oct++) {
                valueMap2D(x * f, y * f, f / np.spread.X, f / np.spread.Y, seed + np.seed + oct, permutations);
                updateResults(g, persist_buf, persistence_map, bufsize);
                f *= np.lacunarity;
                g *= np.persist;
            }

            if (Math.abs(np.offset) > 0.00001f || Math.abs(np.scale - 1.0f) > 0.00001f) {
                for (int i = 0; i < bufsize; i++) {
                    result[i] = result[i] * np.scale + np.offset;
                }
            }

            return result;
        }

        float[] noiseMap3D(float x, float y, float z, float[] persistence_map, int[] permutations) {
            float f = 1.0f, g = 1.0f;
            int bufsize = sx * sy * sz;

            x /= np.spread.X;
            y /= np.spread.Y;
            z /= np.spread.Z;

            for (int i = 0; i < bufsize; i++) {
                result[i] = 0.0f;
            }

            if (persistence_map != null) {
                if (persist_buf == null) {
                    persist_buf = new float[bufsize];
                }
                for (int i = 0; i < bufsize; i++) {
                    persist_buf[i] = 1.0f;
                }
            }

            for (int oct = 0; oct < np.octaves; oct++) {
                valueMap3D(x * f, y * f, z * f, f / np.spread.X, f / np.spread.Y, f / np.spread.Z, seed + np.seed + oct, permutations);
                updateResults(g, persist_buf, persistence_map, bufsize);
                f *= np.lacunarity;
                g *= np.persist;
            }

            if (Math.abs(np.offset) > 0.00001f || Math.abs(np.scale - 1.0f) > 0.00001f) {
                for (int i = 0; i < bufsize; i++) {
                    result[i] = result[i] * np.scale + np.offset;
                }
            }

            return result;
        }

        void valueMap2D(float x, float y, float step_x, float step_y, int seed, int[] permutations) {
            boolean eased = (np.flags & 2) != 0; // NOISE_FLAG_EASED
            int x0 = (int) Math.floor(x);
            int y0 = (int) Math.floor(y);
            float u = x - x0;
            float v = y - y0;
            float orig_u = u;

            int nlx = (int) Math.ceil(u + sx * step_x) + 3;
            int nly = (int) Math.ceil(v + sy * step_y) + 3;
            int index = 0;
            for (int j = 0; j < nly; j++) {
                for (int i = 0; i < nlx; i++) {
                    noise_buf[index++] = noise2d(x0 + i, y0 + j, seed, permutations);
                }
            }

            index = 0;
            int noisey = 0;
            for (int j = 0; j < sy; j++) {
                float v00 = noise_buf[noisey * nlx];
                float v10 = noise_buf[noisey * nlx + 1];
                float v01 = noise_buf[(noisey + 1) * nlx];
                float v11 = noise_buf[(noisey + 1) * nlx + 1];

                u = orig_u;
                int noisex = 0;
                for (int i = 0; i < sx; i++) {
                    value_buf[index++] = biLinearInterpolation(v00, v10, v01, v11, u, v, eased);
                    u += step_x;
                    if (u >= 1.0f) {
                        u -= 1.0f;
                        noisex++;
                        v00 = v10;
                        v01 = v11;
                        v10 = noise_buf[noisey * nlx + noisex + 1];
                        v11 = noise_buf[(noisey + 1) * nlx + noisex + 1];
                    }
                }
                v += step_y;
                if (v >= 1.0f) {
                    v -= 1.0f;
                    noisey++;
                }
            }
        }

        void valueMap3D(float x, float y, float z, float step_x, float step_y, float step_z, int seed, int[] permutations) {
            boolean eased = (np.flags & 2) != 0; // NOISE_FLAG_EASED
            int x0 = (int) Math.floor(x);
            int y0 = (int) Math.floor(y);
            int z0 = (int) Math.floor(z);
            float u = x - x0;
            float v = y - y0;
            float w = z - z0;
            float orig_u = u;
            float orig_v = v;

            int nlx = (int) Math.ceil(u + sx * step_x) + 3;
            int nly = (int) Math.ceil(v + sy * step_y) + 3;
            int nlz = (int) Math.ceil(w + sz * step_z) + 3;

            if (nlx * nly * nlz > noise_buf.length) {
                noise_buf = new float[nlx * nly * nlz];
            }

            int index = 0;
            for (int k = 0; k < nlz; k++) {
                for (int j = 0; j < nly; j++) {
                    for (int i = 0; i < nlx; i++) {
                        noise_buf[index++] = noise3d(x0 + i, y0 + j, z0 + k, seed, permutations);
                    }
                }
            }

            index = 0;
            int noisey = 0;
            int noisez = 0;
            for (int k = 0; k < sz; k++) {
                v = orig_v;
                noisey = 0;
                for (int j = 0; j < sy; j++) {
                    if (noisey >= nly - 1 || noisez >= nlz - 1) {
                        continue;
                    }
                    float v000 = noise_buf[(noisez * nly + noisey) * nlx];
                    float v100 = noise_buf[(noisez * nly + noisey) * nlx + 1];
                    float v010 = noise_buf[(noisez * nly + noisey + 1) * nlx];
                    float v110 = noise_buf[(noisez * nly + noisey + 1) * nlx + 1];
                    float v001 = noise_buf[((noisez + 1) * nly + noisey) * nlx];
                    float v101 = noise_buf[((noisez + 1) * nly + noisey) * nlx + 1];
                    float v011 = noise_buf[((noisez + 1) * nly + noisey + 1) * nlx];
                    float v111 = noise_buf[((noisez + 1) * nly + noisey + 1) * nlx + 1];

                    u = orig_u;
                    int noisex = 0;
                    for (int i = 0; i < sx; i++) {
                        if (index >= value_buf.length) {
                            break;
                        }
                        value_buf[index++] = triLinearInterpolation(v000, v100, v010, v110, v001, v101, v011, v111, u, v, w, eased);
                        u += step_x;
                        if (u >= 1.0f) {
                            u -= 1.0f;
                            noisex++;
                            if (noisex + 1 < nlx && noisey < nly && noisez < nlz) {
                                v000 = v100;
                                v010 = v110;
                                v001 = v101;
                                v011 = v111;
                                v100 = noise_buf[(noisez * nly + noisey) * nlx + noisex + 1];
                                v110 = noise_buf[(noisez * nly + noisey + 1) * nlx + noisex + 1];
                                v101 = noise_buf[((noisez + 1) * nly + noisey) * nlx + noisex + 1];
                                v111 = noise_buf[((noisez + 1) * nly + noisey + 1) * nlx + noisex + 1];
                            }
                        }
                    }
                    v += step_y;
                    if (v >= 1.0f) {
                        v -= 1.0f;
                        noisey++;
                    }
                }
                w += step_z;
                if (w >= 1.0f) {
                    w -= 1.0f;
                    noisez++;
                }
            }
        }

        void updateResults(float g, float[] gmap, float[] persistence_map, int bufsize) {
            boolean absvalue = (np.flags & 4) != 0; // NOISE_FLAG_ABSVALUE
            if (absvalue) {
                if (persistence_map != null) {
                    for (int i = 0; i < bufsize; i++) {
                        result[i] += gmap[i] * Math.abs(value_buf[i]);
                        gmap[i] *= persistence_map[i];
                    }
                } else {
                    for (int i = 0; i < bufsize; i++) {
                        result[i] += g * Math.abs(value_buf[i]);
                    }
                }
            } else {
                if (persistence_map != null) {
                    for (int i = 0; i < bufsize; i++) {
                        result[i] += gmap[i] * value_buf[i];
                        gmap[i] *= persistence_map[i];
                    }
                } else {
                    for (int i = 0; i < bufsize; i++) {
                        result[i] += g * value_buf[i];
                    }
                }
            }
        }
    }

    public NoiseGeneratorPerlin() {
        this(new Random());
    }

    public NoiseGeneratorPerlin(Random var1) {
        this.permutations = new int[512];
        this.xCoord = var1.nextDouble() * 256.0D;
        this.yCoord = var1.nextDouble() * 256.0D;
        this.zCoord = var1.nextDouble() * 256.0D;

        for (int var2 = 0; var2 < 256; this.permutations[var2] = var2++) {
        }

        for (int var2 = 0; var2 < 256; ++var2) {
            int var3 = var1.nextInt(256 - var2) + var2;
            int var4 = this.permutations[var2];
            this.permutations[var2] = this.permutations[var3];
            this.permutations[var3] = var4;
            this.permutations[var2 + 256] = this.permutations[var2];
        }

        NoiseParams np = new NoiseParams();
        np.seed = var1.nextInt();
        this.noise = new Noise(np, np.seed, 1, 1, 1);
    }

    public double generateNoise(double var1, double var3, double var5) {
        float value = noise.noise3d_value((float) (var1 + xCoord), (float) (var3 + yCoord), (float) (var5 + zCoord),
                noise.seed, (noise.np.flags & 2) != 0, permutations);
        return value;
    }

    public double func_801_a(double var1, double var3) {
        float value = noise.noise2d_value((float) (var1 + xCoord), (float) (var3 + yCoord),
                noise.seed, (noise.np.flags & 2) != 0, permutations);
        return value;
    }

    public void func_805_a(double[] var1, double var2, double var4, double var6, int var8, int var9, int var10,
            double var11, double var13, double var15, double var17) {
        noise = new Noise(noise.np, noise.seed, var8, var9, var10);

        if (var9 == 1) {
            float[] result = noise.noiseMap2D((float) var2, (float) var6, null, permutations);
            for (int i = 0; i < var8 * var10; i++) {
                var1[i] += result[i] / var17;
            }
        } else {
            float[] result = noise.noiseMap3D((float) var2, (float) var4, (float) var6, null, permutations);
            for (int i = 0; i < var8 * var9 * var10; i++) {
                var1[i] += result[i] / var17;
            }
        }
    }

    public final double lerp(double var1, double var3, double var5) {
        return var3 + var1 * (var5 - var3);
    }

    public final double func_4110_a(int var1, double var2, double var4) {
        int var6 = var1 & 15;
        double var7 = (double) (1 - ((var6 & 8) >> 3)) * var2;
        double var9 = var6 < 4 ? 0.0D : (var6 != 12 && var6 != 14 ? var4 : var2);
        return ((var6 & 1) == 0 ? var7 : -var7) + ((var6 & 2) == 0 ? var9 : -var9);
    }

    public final double grad(int var1, double var2, double var4, double var6) {
        int var8 = var1 & 15;
        double var9 = var8 < 8 ? var2 : var4;
        double var11 = var8 < 4 ? var4 : (var8 != 12 && var8 != 14 ? var6 : var2);
        return ((var8 & 1) == 0 ? var9 : -var9) + ((var8 & 2) == 0 ? var11 : -var11);
    }
}