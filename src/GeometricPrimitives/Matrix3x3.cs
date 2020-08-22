using System;
using MGSharp.Core.GeometricPrimitives;

namespace MGSharp.Core.LinearAlgebra
{
    public class Matrix3x3
    {
        // Constants
        const double pi = 3.14159265358979323846;
        const int outprecision = 17;
        const int outwidth = 26;

        // dlambda_limit, below which two lambdas are relatively equal
        double dlambda_limit = 1.0E-3;
        double iszero_limit = 1.0E-20;

        // Some globals to record global statistics
        int n_all_lambdas_equal = 0, n_two_lambdas_equal = 0,
            n_all_lambdas_different = 0;

        private double[,] mat;
        public Vector eigenValues;
        public Vector eulerAngles;
        public Matrix3x3 R;
        public Matrix3x3 D;

        //constructor
        public Matrix3x3()
        {
            //create identity matrix
            int row, col;
            mat = new double[3, 3];
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    if (row == col)
                        mat[row, col] = 0.0;
                    else
                        mat[row, col] = 0.0;
                }
            }
        }

        public Matrix3x3(double v11, double v12, double v13, 
                        double v21, double v22, double v23,
                            double v31, double v32, double v33)
        {
            mat = new double[3, 3];
            
            mat[0, 0] = v11; mat[0, 1] = v12; mat[0, 2] = v13;
            mat[1, 0] = v21; mat[1, 1] = v22; mat[1, 2] = v23;
            mat[2, 0] = v31; mat[2, 1] = v32; mat[2, 2] = v33;

            eigenValues = new Vector();
            eulerAngles = new Vector();
            R = new Matrix3x3();
            D = new Matrix3x3();
        }

        public Matrix3x3(Vector u, double theta)
        {
            mat = new double[3, 3];

            double cos = Math.Cos(theta);
            double sin = Math.Sin(theta);

            mat[0, 0] = cos + u.x * u.x * (1 - cos);
            mat[0, 1] = - u.z * sin + u.x * u.y * (1 - cos);
            mat[0, 2] = u.y * sin + u.x * u.z * (1 - cos);
            mat[1, 0] = u.z * sin + u.x * u.y * (1 - cos);
            mat[1, 1] = cos + u.y * u.y * (1 - cos);
            mat[1, 2] = - u.x * sin + u.y * u.z * (1 - cos);
            mat[2, 0] = - u.y * sin + u.x * u.z * (1 - cos);
            mat[2, 1] = u.x * sin + u.y * u.z * (1 - cos);
            mat[2, 2] = cos + u.z * u.z * (1 - cos);
        }

        public void SetZero()
        {
            int row, col;
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    if (row == col)
                        mat[row, col] = 0.0;
                    else
                        mat[row, col] = 0.0;
                }
            }
        }

        public object Clone()
        {
            Matrix3x3 m = new Matrix3x3();
            int row, col;
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    m[row, col] = mat[row, col];
                }
            }
            return (object)m;
        }

        public override bool Equals(object b)
        {
            Matrix3x3 m = (Matrix3x3)b;
            int row, col;
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    if (m[row, col] != mat[row, col]) return false;
                }
            }
            return true;
        }

        //index operator, to retrieve and set values 
        //  in the matrix
        public double this[int row, int col]
        {
            get { return mat[row, col]; }
            set { mat[row, col] = value; }
        }

        //arithmetic operators
        public static Matrix3x3 operator +(Matrix3x3 a, Matrix3x3 b)
        {
            int row, col;
            Matrix3x3 r = new Matrix3x3();
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    r[row, col] = a[row, col] + b[row, col];
                }
            }
            return r;
        }

        public static Matrix3x3 operator -(Matrix3x3 a, Matrix3x3 b)
        {
            int row, col;
            Matrix3x3 r = new Matrix3x3();
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    r[row, col] = a[row, col] - b[row, col];
                }
            }
            return r;
        }

        public static Matrix3x3 operator *(Matrix3x3 a, Matrix3x3 b)
        {
            int row, col, i;
            Matrix3x3 r = new Matrix3x3();
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    r[row, col] = 0;
                    for (i = 0; i < 3; i++)
                    {
                        r[row, col] += a[row, i] * b[i, col];
                    }
                }
            }
            return r;
        }

        public static Vector operator *(Matrix3x3 a, Vector v)
        {
            Vector r = new Vector();
            r.x = v.x * a[0, 0] + v.y * a[0, 1] + v.z * a[0, 2];// + v.w * a[0, 3];
            r.y = v.x * a[1, 0] + v.y * a[1, 1] + v.z * a[1, 2];// + v.w * a[1, 3];
            r.z = v.x * a[2, 0] + v.y * a[2, 1] + v.z * a[2, 2];// + v.w * a[2, 3];
            //r.w = v.x * a[3, 0] + v.y * a[3, 1] + v.z * a[3, 2] + v.w * a[3, 3];
            r.w = 0;

            return r;
        }

        public Matrix3x3 Transpose()
        {
            int row, col;
            Matrix3x3 m = new Matrix3x3();
            for (row = 0; row < 3; row++)
            {
                for (col = 0; col < 3; col++)
                {
                    m[row, col] = mat[col, row];
                }
            }
            return m;
        }

        public static Matrix3x3 Translate(Vector t)
        {
            Matrix3x3 m = new Matrix3x3();
            m[0, 3] = t.x;
            m[1, 3] = t.y;
            m[2, 3] = t.z;
            return m;
        }

        public static Matrix3x3 Scale(double s)
        {
            Matrix3x3 m = new Matrix3x3();
            m[0, 0] = s;
            m[1, 1] = s;
            m[2, 2] = s;
            return m;
        }

        public static Matrix3x3 Scale(double sx, double sy, double sz)
        {
            Matrix3x3 m = new Matrix3x3();
            m[0, 0] = sx;
            m[1, 1] = sy;
            m[2, 2] = sz;
            return m;
        }

        public static Matrix3x3 InvertDiagonalMatrix(Matrix3x3 diag)
        {
            Matrix3x3 diagInverse = new Matrix3x3();

            diagInverse.mat[0, 0] = 1 / diag[0, 0];
            diagInverse.mat[1, 1] = 1 / diag[1, 1];
            diagInverse.mat[2, 2] = 1 / diag[2, 2];

            return diagInverse;
        }

        //general rotation around the vector axis
        public static Matrix3x3 Rotate(Vector axis, double angle)
        {
            Matrix3x3 r = new Matrix3x3();
            double c = Math.Cos(angle);
            double s = Math.Sin(angle);

            axis.normalize();
            r[0, 0] = (1 - c) * axis.x * axis.x + c;
            r[1, 0] = (1 - c) * axis.x * axis.y + s * axis.z;
            r[2, 0] = (1 - c) * axis.x * axis.z - s * axis.y;
            r[3, 0] = 0;

            r[0, 1] = (1 - c) * axis.x * axis.y - s * axis.z;
            r[1, 1] = (1 - c) * axis.y * axis.y + c;
            r[2, 1] = (1 - c) * axis.y * axis.z + s * axis.x;
            r[3, 1] = 0;

            r[0, 2] = (1 - c) * axis.x * axis.z + s * axis.y;
            r[1, 2] = (1 - c) * axis.y * axis.z - s * axis.x;
            r[2, 2] = (1 - c) * axis.z * axis.z + c;
            r[3, 2] = 0;

            r[0, 3] = 0;
            r[1, 3] = 0;
            r[2, 3] = 0;
            r[3, 3] = 1;

            return r;
        }

        public static Matrix3x3 Projection(double n, double f,
            double t, double b,
            double l, double r)
        {
            Matrix3x3 m = new Matrix3x3();

            m[0, 0] = 2 * n / (r - l);
            m[1, 0] = 0;
            m[2, 0] = 0;
            m[3, 0] = 0;

            m[0, 1] = 0;
            m[1, 1] = 2 * n / (t - b);
            m[2, 1] = 0;
            m[3, 1] = 0;

            m[0, 2] = 0;
            m[1, 2] = 2 * n / (t - b);
            m[2, 2] = -(f + n) / (f - n);
            m[3, 2] = -1;

            m[0, 3] = (r + l) / (r - l);
            m[1, 3] = (t + b) / (t - b);
            m[2, 3] = -2 * f * n / (f - n);
            m[3, 3] = 0;

            return m;
        }

        //provide a method to display the matrix	
        public override string ToString()
        {
            int i;
            String s = "\n";
            for (i = 0; i < 3; i++)
            {
                //s += String.Format("[{0},{1},{2},{3}]\n",
                //    mat[i, 0], mat[i, 1], mat[i, 2], mat[i, 3]);

                s += String.Format("[{0},{1},{2}]\n",
                    mat[i, 0], mat[i, 1], mat[i, 2]);
            }
            return s;
        }

        double sqr(double val)
        {
            return val * val;
        }

        double angle(double x, double y)
        {
            if (x == 0.0)
                return (y == 0.0 ? 0.0 : 0.5 * pi * Math.Sign(y));
            return (x < 0.0 ? Math.Atan(y / x) + pi * Math.Sign(y)
                            : Math.Atan(y / x));
        }

        // The functions trunc_sqrt and trunc_acos prevent a domain error
        // because of rounding off errors:

        double trunc_sqrt(double x)
        {
            return (x <= 0.0 ? 0.0 : Math.Sqrt(x));
        }

        double trunc_acos(double x)
        {
            if (x >= 1.0)
                return 0.0;
            if (x <= -1.0)
                return pi;
            return Math.Acos(x);
        }

        // Solve the lambdas from the matrix:
        void solve_lambdas()
        {
            double p, q, b, t, delta;
            b = mat[0,0] + mat[1,1] + mat[2,2];
            t = sqr(mat[0,1]) + sqr(mat[0,2]) + sqr(mat[1,2]);
            p = 0.5 * (sqr(mat[0,0] - mat[1,1]) + sqr(mat[0,0] - mat[2,2])
                + sqr(mat[1,1] - mat[2,2]));
            p += 3.0 * t;
            q = 18.0 * (mat[0,0] * mat[1,1] * mat[2,2] + 3.0 * mat[0,1] * mat[0,2] * mat[1,2]);
            q += 2.0 * (mat[0,0] * sqr(mat[0,0]) + mat[1,1] * sqr(mat[1,1]) +
                         mat[2,2] * sqr(mat[2,2]));

            q += 9.0 * b * t;
            q -= 3.0 * (mat[0,0] + mat[1,1]) * (mat[0,0] + mat[2,2]) *
                       (mat[1,1] + mat[2,2]);
            q -= 27.0 * (mat[0,0] * sqr(mat[1,2]) + mat[1,1] * sqr(mat[0,2]) +
                          mat[2,2] * sqr(mat[0,1]));
            if (p < iszero_limit)
                eigenValues.x = eigenValues.y = eigenValues.z = b / 3.0;
            else
            {
                delta = trunc_acos(0.5 * q / Math.Sqrt(p * sqr(p)));
                p = 2.0 * Math.Sqrt(p);
                // Changing the order in result yields different angles but identical matrix
                eigenValues.x = (b + p * Math.Cos(delta / 3.0)) / 3.0;
                eigenValues.y = (b + p * Math.Cos((delta + 2.0 * pi) / 3.0)) / 3.0;
                eigenValues.z = (b + p * Math.Cos((delta - 2.0 * pi) / 3.0)) / 3.0;
            };
        }

        // Determine which type of solution is needed:
        //  0: all lambdas equal
        //  1: two lambdas equal
        //  2: all lambdas different

        int solve_type(double[] lambdas)
        {
            int i1 = 0, i2 = 0, isum = 0;
            double t, lambdasum = 0.0;

            for (int i = 0; i < 3; i++)
                lambdasum += sqr(lambdas[i]);

            lambdasum = Math.Sqrt(lambdasum);
            for (int i = 0; i < 2; i++)
                for (int j = i + 1; j < 3; j++)
                {
                    t = Math.Abs(lambdas[i] - lambdas[j]);
                    if (lambdasum > iszero_limit)
                        t /= lambdasum;
                    if (t < dlambda_limit)
                    {
                        isum++;
                        i1 = i;
                        i2 = j;
                    };
                };
            if (isum == 0)
                return 2;
            if (isum >= 2)
                return 0;
            t = 0.5 * (lambdas[i1] + lambdas[i2]);
            lambdas[2] = lambdas[3 - i1 - i2];
            lambdas[0] = lambdas[1] = t;

            eigenValues.x = lambdas[0];
            eigenValues.y = lambdas[1];
            eigenValues.z = lambdas[2];

            return 1;
        }

        void set_diagonal(Vector arg)
        {
            D[0, 0] = arg.x;
            D[1, 1] = arg.y;
            D[2, 2] = arg.z;
        }

        // Solve the angles from the matrix and the solved lambdas:
        //  solve_angles_0: all lambdas equal
        //  solve_angles_1: two lambdas equal
        //  solve_angles_2: all lambdas different
        public Matrix3x3 sym_eigen()//Vector lambdas, Vector angles)
        {
            solve_lambdas();
            switch(solve_type(eigenValues.ToArray()) )
            { case 0:
                solve_angles_0(eulerAngles.ToArray(), eigenValues.ToArray());
                n_all_lambdas_equal++;
                break;
            case 1:
                solve_angles_1(eulerAngles.ToArray(), eigenValues.ToArray());
                n_two_lambdas_equal++;
                break;
            case 2:
                solve_angles_2(eulerAngles.ToArray(), eigenValues.ToArray());
                n_all_lambdas_different++;
                break;
            };

            set_diagonal(eigenValues);

            return this;
        }

        public void set_rotation(int axis, double phi )
        {
            if (axis == 0 || axis > 3 )
                throw new Exception( "set_rotation: axis number out of bounds." );

            SetZero();
            
            switch (axis )
            {
                case 1:
                    mat[0,0] = 1.0;
                    mat[1,1] = mat[2,2] = Math.Cos(phi);
                    mat[1,2] = mat[2,1] = Math.Sin(phi);
                    mat[1,2] = -mat[1,2];
                    break;
                case 2:
                    mat[1,1] = 1.0;
                    mat[0,0] = mat[2,2] = Math.Cos(phi);
                    mat[0,2] = mat[2,0] = Math.Sin(phi);
                    mat[2,0] = -mat[2,0];
                    break;
                case 3:
                    mat[2,2] = 1.0;
                    mat[0,0] = mat[1,1] = Math.Cos(phi);
                    mat[0,1] = mat[1,0] = Math.Sin(phi);
                    mat[0,1] = -mat[0,1];
                    break;
            }
        }

        // Generate the symmetric matrix A from lambdas and angles:
        public void compute_matrix()//Vector lambdas, Vector angles )
        {
            Matrix3x3 rot = new Matrix3x3(), temp = new Matrix3x3();

            //Console.WriteLine("Orthogonal matrix (Empty):");
            //Console.WriteLine(R);

            R.set_rotation( 3, eulerAngles.z );

            //Console.WriteLine("Orthogonal matrix (1st rotation):");
            //Console.WriteLine(R);

            rot.set_rotation( 2, eulerAngles.y );
            rot *= R;
            R = (Matrix3x3)rot.Clone();

            //Console.WriteLine("Orthogonal matrix (2nd rotation):");
            //Console.WriteLine(rot);
            //Console.WriteLine(R);

            rot.set_rotation( 1, eulerAngles.x );
            //R *= rot;
            rot *= R;
            R = (Matrix3x3)rot.Clone();

            //Console.WriteLine("Orthogonal matrix (3rd rotation):");
            //Console.WriteLine(R);
            //Console.WriteLine(rot);

            //
            //temp = R;
            //R.Transpose();

            //rot.SetZero();
            //rot[0, 0] = eigenValues.x;
            //rot[1, 1] = eigenValues.y;
            //rot[2, 2] = eigenValues.z;
            //rot.set_diagonal(eigenValues );


            //R *= rot;
            //R *= temp;
        }

        void solve_angles_0(double[] res, double[] lambdas )
        {
            res[0] = 0.0;
            res[1] = 0.0;
            res[2] = 0.0;

            eulerAngles = new Vector();
        }

        void solve_angles_1(double[] res, double[] lambdas )
        {
            double phi1a, phi1b, phi2, absdif, delta = 1.0E10;
            double g12, g21, t1, t2;
            t1 = lambdas[0] - lambdas[2];
            t2 = mat[0,0] - lambdas[2];
            phi2 = trunc_acos(trunc_sqrt(t2 / t1)); // + pi for symmetry
            g12 = 0.5 * t1 * Math.Sin(2.0 * phi2);
            g21 = t2;
            t1 = angle(mat[1,1] - mat[2,2], -2.0 * mat[1,2]);
            t2 = Math.Sin(t1);
            t1 = Math.Cos(t1);
            phi1b = 0.5 * angle(g21 * t1, -g21 * t2);
            t1 = angle(mat[0,1], -1.0 * mat[0,2]);
            t2 = Math.Sin(t1);
            t1 = Math.Cos(t1);
            bool big = sqr(mat[1,1] - mat[2,2]) + sqr(2.0 * mat[1,2])
                      > sqr(mat[0,1]) + sqr(mat[0,2]);
            for (int i = 0; i < 2; i++)
            {
                phi1a = angle(g12 * t2, g12 * t1);
                absdif = Math.Abs(phi1a - phi1b);
                if (absdif < delta)
                {
                    delta = absdif;
                    res[0] = big ? phi1b : phi1a;
                    res[1] = phi2;
                };
                phi2 = -phi2;
                g12 = -g12;
            };
            res[2] = 0.0;

            eulerAngles.x = res[0];
            eulerAngles.y = res[1];
            eulerAngles.z = res[2];
        }

        void solve_angles_2(double[] res, double[] lambdas )
        {
            double phi1a, phi1b, phi2, phi3, v, w, absdif, delta = 1.0E10;
            double g11, g12, g21, g22, t1, t2, t3, t4;
            t1 = lambdas[0] - lambdas[1];
            t2 = lambdas[1] - lambdas[2];
            t3 = lambdas[2] - lambdas[0];
            t4 = mat[0,0] - lambdas[2];
            v = sqr(mat[0,1]) + sqr(mat[0,2]);
            v += t4 * (mat[0,0] + t3 - lambdas[1]);
            v /= t2 * t3;
            if (Math.Abs(v) < iszero_limit) w = 1.0;
            else w = (t4 - t2 * v) / (t1 * v);
            phi2 = trunc_acos(trunc_sqrt(v)); // + pi for symmetry
            phi3 = trunc_acos(trunc_sqrt(w)); // + pi for symmetry
            g11 = 0.5 * t1 * Math.Cos(phi2) * Math.Sin(2.0 * phi3);
            g12 = 0.5 * (t1 * w + t2) * Math.Sin(2.0 * phi2);
            g21 = t1 * (1.0 + (v - 2.0) * w) + t2 * v;
            g22 = t1 * Math.Sin(phi2) * Math.Sin(2.0 * phi3);
            t1 = angle(mat[0,1], -1.0 * mat[0,2]);
            t3 = angle(mat[1,1] - mat[2,2], -2.0 * mat[1,2]);
            t2 = Math.Sin(t1);
            t1 = Math.Cos(t1);
            t4 = Math.Sin(t3);
            t3 = Math.Cos(t3);
            bool big = sqr(mat[1,1] - mat[2,2]) + sqr(2.0 * mat[1,2])
                      > sqr(mat[0,1]) + sqr(mat[0,2]);
            for (int i = 0; i < 4; i++)
            {
                phi1a = angle(g11 * t1 + g12 * t2, -g11 * t2 + g12 * t1);
                phi1b = 0.5 * angle(g21 * t3 + g22 * t4, -g21 * t4 + g22 * t3);
                absdif = Math.Abs(phi1a - phi1b);
                if (absdif < delta)
                {
                    delta = absdif;
                    res[0] = big ? phi1b : phi1a;
                    res[1] = phi2;
                    res[2] = phi3;
                };
                phi3 = -phi3;
                g11 = -g11;
                g22 = -g22;
                if (i == 1)
                {
                    phi2 = -phi2;
                    g12 = -g12;
                    g22 = -g22;
                }
            }

            eulerAngles.x = res[0];
            eulerAngles.y = res[1];
            eulerAngles.z = res[2];
        }

        public static void Main(string[] args)
        {
            Matrix3x3 input_matrix = new Matrix3x3(1, 2, 3, 2, 4, 5, 3, 5, 6);
            
            input_matrix.sym_eigen();
            Console.WriteLine("Original matrix:");
            Console.WriteLine(input_matrix);

            Console.WriteLine("Diagonal matrix:");
            Console.WriteLine(input_matrix.D);

            input_matrix.compute_matrix();
            Console.WriteLine("Orthogonal matrix:");
            Console.WriteLine(input_matrix.R);

            Console.WriteLine("Solved matrix:");
            //Console.WriteLine(input_matrix.R.Transpose() * input_matrix.R);
            Console.WriteLine(input_matrix.R * input_matrix.D * input_matrix.R.Transpose());

            //input_matrix.eulerAngles *= 180 / pi;

            //Console.WriteLine("Anti-Clockwise Euler Angles (degrees): ");
            //Console.WriteLine(input_matrix.eulerAngles);

            Console.ReadLine();
        }
    };
}