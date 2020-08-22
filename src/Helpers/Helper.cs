using System;
using System.Collections;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;

namespace MGSharp.Core.Helpers
{
    public class Helper
    {
        public static Vector[] NCubedInN;
        
        public static void SetNCubedInN(int size)
        {
            NCubedInN = new Vector[size];
            int dim = (int) Math.Round(Math.Pow(size, 1 / 3f));
            int i = 0;
            for (int z = 0; z < dim; z++)
            {
                for (int y = 0; y < dim; y++)
                {
                    for (int x = 0; x < dim; x++)
                    {
                        NCubedInN[i] = new Vector(x, z, y);
                        //Console.WriteLine(i + ": " + NCubedInN[i]);
                        i++;
                    }
                }
            }
        }

        public static void SetNCubedInN(Vector dim)
        {
            NCubedInN = new Vector[(int) (dim.x*dim.y*dim.z)];
            int dimX = (int)dim.x;
            int dimY = (int)dim.y;
            int dimZ = (int)dim.z;

            int i = 0;

            for (int y = 0; y < dimY; y++)
            {
                for (int z = 0; z < dimZ; z++)
                {
                    for (int x = 0; x < dimX; x++)
                    {
                        NCubedInN[i] = new Vector(x, y, z);
                        //Console.WriteLine(i + ": " + NCubedInN[i]);
                        i++;
                    }
                }
            }
        }

        public static List<Edge> Distinct(List<Edge> edges)
        {
            List<Edge> newEdges = new List<Edge>();
            Console.WriteLine("Before: " + edges.Count);
            for(int i=0; i<edges.Count; i++)
            {
                newEdges.Add(edges[i]);
                for(int j=i+1; j<edges.Count; j++)
                {
                    if (edges[i].ends[0].pos == edges[j].ends[0].pos &&
                        edges[i].ends[1].pos == edges[j].ends[1].pos &&
                        edges[i].ends[0].cellId == edges[j].ends[0].cellId &&
                        edges[i].ends[1].cellId == edges[j].ends[1].cellId)
                    {
                        newEdges.Remove(edges[i]);
                    }
                }
            }
            Console.WriteLine("After: " + newEdges.Count);
            return newEdges;
        }

        public static int[] OrderEggCylinderCells(int maxSize, int popSize)
        {
            List<Vector> NCubedInN1 = new List<Vector>(maxSize);

            List<int> boundary = new List<int>();
            int dim = (int)Math.Floor(Math.Pow(maxSize, 1f / 3f));
            int dim2 = dim * dim;

            for (int i=0; i<popSize; i++)
            {
                int y = i / dim2;
                int xz = i % dim2;
                int x = xz / dim;
                int z = xz % dim;

                if(x==0 || x == dim || y == 0 || y == dim || z == 0 || z == dim)
                {
                    boundary.Add(i);
                    NCubedInN1.Add(new Vector(x, y, z));
                }
            }

            for (int i = 0; i < popSize; i++)
            {
                int y = i / dim2;
                int xz = i % dim2;
                int x = xz / dim;
                int z = xz % dim;

                if (!(x == 0 || x == dim || y == 0 || y == dim || z == 0 || z == dim))
                {
                    boundary.Add(i);
                    NCubedInN1.Add(new Vector(x, y, z));
                }
            }
            
            NCubedInN = NCubedInN1.ToArray();
            //Console.Write(NCubedInN[171]);
            return boundary.ToArray();
        }

        public static List<int> GetMooreNeighbourhood3D(int cellIndex, int maxSize, int size)
        {
            int dim = (int)Math.Floor(Math.Pow(maxSize, 1f / 3f));
            int dim2 = dim * dim;

            int z = cellIndex / dim2;
            int xy = cellIndex % dim2;
            int x = xy / dim;
            int y = xy % dim;

            int[] neighbours = {
	            (cellIndex - dim2 - dim - 1),
	            (cellIndex - dim2 - dim),
	            (cellIndex - dim2 - dim + 1),
	            (cellIndex - dim2 - 1),
	            (cellIndex - dim2),
	            (cellIndex - dim2 + 1),
	            (cellIndex - dim2 + dim - 1),
	            (cellIndex - dim2 + dim),
	            (cellIndex - dim2 + dim + 1),
                (cellIndex - dim - 1) ,
                (cellIndex - dim),
                (cellIndex - dim + 1),
                (cellIndex - 1),
                (cellIndex + 1),
                (cellIndex + dim - 1),
                (cellIndex + dim),
                (cellIndex + dim + 1),
                (cellIndex + dim2 - dim - 1),
                (cellIndex + dim2 - dim),
                (cellIndex + dim2 - dim + 1),
                (cellIndex + dim2 - 1),
                (cellIndex + dim2),
                (cellIndex + dim2 + 1),
                (cellIndex + dim2 + dim - 1),
                (cellIndex + dim2 + dim),
                (cellIndex + dim2 + dim + 1),
            };

            List<int> mooreNeighbours = new List<int>();

            for (int j = 0; j < 9; j++)
            {
                int i = j;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z - 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 6 && nz == z - 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 6 && nz == z - 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }
            for (int j = 9; j < 17; j++)
            {
                int i = j - 9;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 5 && nz == z && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 5 && nz == z && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }
            for (int j = 17; j < 26; j++)
            {
                int i = j - 17;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z + 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 6 && nz == z + 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 6 && nz == z + 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }

            return mooreNeighbours;
        }

        public static List<int> GetMooreNeighbourhood3D(int cellIndex, Vector maxSize, int size)
        {
            int dimX = (int) maxSize.x;
            int dimY = (int)maxSize.y;
            int dimZ = (int)maxSize.z;

            int dim2 = dimX * dimZ;

            int z = cellIndex / dim2;
            int xy = cellIndex % dim2;
            int x = xy / dimZ;
            int y = xy % dimZ;

            int[] neighbours = {
                (cellIndex - dim2 - dimZ - 1),
                (cellIndex - dim2 - dimZ),
                (cellIndex - dim2 - dimZ + 1),
                (cellIndex - dim2 - 1),
                (cellIndex - dim2),
                (cellIndex - dim2 + 1),
                (cellIndex - dim2 + dimZ - 1),
                (cellIndex - dim2 + dimZ),
                (cellIndex - dim2 + dimZ + 1),
                (cellIndex - dimZ - 1) ,
                (cellIndex - dimZ),
                (cellIndex - dimZ + 1),
                (cellIndex - 1),
                (cellIndex + 1),
                (cellIndex + dimZ - 1),
                (cellIndex + dimZ),
                (cellIndex + dimZ + 1),
                (cellIndex + dim2 - dimZ - 1),
                (cellIndex + dim2 - dimZ),
                (cellIndex + dim2 - dimZ + 1),
                (cellIndex + dim2 - 1),
                (cellIndex + dim2),
                (cellIndex + dim2 + 1),
                (cellIndex + dim2 + dimZ - 1),
                (cellIndex + dim2 + dimZ),
                (cellIndex + dim2 + dimZ + 1),
            };

            List<int> mooreNeighbours = new List<int>();

            for (int j = 0; j < 9; j++)
            {
                int i = j;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dimZ;
                int ny = nxy % dimZ;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z - 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 6 && nz == z - 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 6 && nz == z - 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }
            for (int j = 9; j < 17; j++)
            {
                int i = j - 9;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dimZ;
                int ny = nxy % dimZ;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 5 && nz == z && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 5 && nz == z && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }
            for (int j = 17; j < 26; j++)
            {
                int i = j - 17;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dimZ;
                int ny = nxy % dimZ;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 3 && nz == z + 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 3 && i < 6 && nz == z + 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 6 && nz == z + 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }

            return mooreNeighbours;
        }

        public static void SetNCubedInUniqueHexa()
        {
            float radius = 1f;// (float) (2/Math.Sqrt(3));
            NCubedInN = new Vector[7];

            NCubedInN[0] = Vector.zero;
            for (int i = 1; i < 7; i++)
            {
                float angle = (float) (2 * Math.PI * (i) / 6);
                NCubedInN[i] = new Vector(radius * Math.Cos(angle), 0, radius * Math.Sin(angle));
            }
        }

        public static List<int> GetHexagonalNeighbourhood3D(int cellIndex, int maxSize, int size)
        {
            int dim = (int)Math.Floor(Math.Pow(maxSize, 1f / 3f));
            int dim2 = dim * dim;

            int z = cellIndex / dim2;
            int xy = cellIndex % dim2;
            int x = xy / dim;
            int y = xy % dim;
            int[] neighbours = new int[20];
            if (x % 2 == 1)
            {
                neighbours[0] = (cellIndex - dim2 - dim);
                neighbours[1] = (cellIndex - dim2 - dim + 1);
                neighbours[2] = (cellIndex - dim2 - 1);
                neighbours[3] = (cellIndex - dim2);
                neighbours[4] = (cellIndex - dim2 + 1);
                neighbours[5] = (cellIndex - dim2 + dim);
                neighbours[6] = (cellIndex - dim2 + dim + 1);
                neighbours[7] = (cellIndex - dim);
                neighbours[8] = (cellIndex - dim + 1);
                neighbours[9] = (cellIndex - 1);
                neighbours[10] = (cellIndex + 1);
                neighbours[11] = (cellIndex + dim);
                neighbours[12] = (cellIndex + dim + 1);
                neighbours[13] = (cellIndex + dim2 - dim);
                neighbours[14] = (cellIndex + dim2 - dim + 1);
                neighbours[15] = (cellIndex + dim2 - 1);
                neighbours[16] = (cellIndex + dim2);
                neighbours[17] = (cellIndex + dim2 + 1);
                neighbours[18] = (cellIndex + dim2 + dim);
                neighbours[19] = (cellIndex + dim2 + dim + 1);
            }
            if (x % 2 == 0)
            {
                neighbours[0] = (cellIndex - dim2 - dim-1);
                neighbours[1] = (cellIndex - dim2 - dim);
                neighbours[2] = (cellIndex - dim2 - 1);
                neighbours[3] = (cellIndex - dim2);
                neighbours[4] = (cellIndex - dim2 + 1);
                neighbours[5] = (cellIndex - dim2 + dim-1);
                neighbours[6] = (cellIndex - dim2 + dim);
                neighbours[7] = (cellIndex - dim-1);
                neighbours[8] = (cellIndex - dim);
                neighbours[9] = (cellIndex - 1);
                neighbours[10] = (cellIndex + 1);
                neighbours[11] = (cellIndex + dim-1);
                neighbours[12] = (cellIndex + dim);
                neighbours[13] = (cellIndex + dim2 - dim-1);
                neighbours[14] = (cellIndex + dim2 - dim);
                neighbours[15] = (cellIndex + dim2 - 1);
                neighbours[16] = (cellIndex + dim2);
                neighbours[17] = (cellIndex + dim2 + 1);
                neighbours[18] = (cellIndex + dim2 + dim - 1);
                neighbours[19] = (cellIndex + dim2 + dim);
            }

            List<int> mooreNeighbours = new List<int>();

            for (int j = 0; j < 7; j++)
            {
                int i = j;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 2 && nz == z - 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 2 && i < 5 && nz == z - 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 5 && nz == z - 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }

            for (int j = 7; j < 13; j++)
            {
                int i = j - 7;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 2 && nz == z && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 2 && i < 4 && nz == z && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 4 && nz == z && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }

            for (int j = 13; j < 20; j++)
            {
                int i = j - 13;
                int nz = neighbours[j] / dim2;
                int nxy = neighbours[j] % dim2;
                int nx = nxy / dim;
                int ny = nxy % dim;
                if (neighbours[j] >= 0 && neighbours[j] < size)
                {
                    if (i < 2 && nz == z + 1 && nx == x - 1)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 2 && i < 5 && nz == z + 1 && nx == x)
                        mooreNeighbours.Add(neighbours[j]);
                    else if (i >= 5 && nz == z + 1 && nx == x + 1)
                        mooreNeighbours.Add(neighbours[j]);
                }
            }

            return mooreNeighbours;
        }
        
        public static void SetNCubedInHexa(int size)
        {
            NCubedInN = new Vector[size];
            int dim = (int)Math.Pow(size, 1 / 3f);
            int halfDim = dim / 2;
            int i = 0;
            for (int y = 0; y < dim; y++)
            {
                for (int z = 0; z < dim; z++)
                {
                    for (int x = 0; x < dim; x++)
                    {
                        //NCubedInN[i] = new Vector(x + y * 0.5f - y / 2, y, z);
                        NCubedInN[i] = new Vector(x + z * .5f - z / 2, y, z);
                        i++;
                    }
                }
            }
        }

        public static void SetNCubedInSphere(int size)
        {
            NCubedInN = new Vector[size];
            int dim = (int)Math.Pow(size, 1 / 3f);

            int i = 0;
            for (int z = 0; z < dim; z++)
            {
                for (int y = 0; y < dim; y++)
                {
                    for (int x = 0; x < dim; x++)
                    {
                        NCubedInN[i] = new Vector(x + y * 0.5f - y / 2, y, z);
                        i++;
                    }
                }
            }
        }
        
        public static int[] DetermineElongationAxis(float[,] table, int length)
        {
            int[] axis = new int[2];
            float max = table[0, 0];

            for (int i = 0; i < length; ++i)
            {
                int j = 0;
                while (j < i)
                {
                    if (table[i, j] > max)
                    {
                        max = table[i, j];
                        axis[0] = i;
                        axis[1] = j;
                    }
                    j++;
                }
            }
            return axis;
        }

        public static int RoundToInt(double x)
        {
            if(x-Math.Floor(x) < 0.5)
            {
                return (int)Math.Floor(x);
            }
            else
            {
                return (int)Math.Ceiling(x);
            }
        }
        
        static int seed = 0;
        public static double GaussianNumber(double mean, double stdDev)
        {
            Random rand = new Random(seed); //reuse this if you are generating many
            seed++;

            double u1 = 1.0 - rand.NextDouble(); //uniform(0,1] random doubles
            double u2 = 1.0 - rand.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            double randNormal =
                         mean + stdDev * randStdNormal; //random normal(mean,stdDev^2)

            return randNormal;
        }
        
        public static void Sort(ref int[] inputArray, ref int[] indices, int n)
        {
            for (int i = 0; i < n - 1; i++)
            {
                for (int j = i + 1; j > 0; j--)
                {
                    if (inputArray[j - 1] > inputArray[j])
                    {
                        int temp = inputArray[j - 1];
                        inputArray[j - 1] = inputArray[j];
                        inputArray[j] = temp;

                        temp = indices[j - 1];
                        indices[j - 1] = indices[j];
                        indices[j] = temp;
                    }
                }
            }
        }

        public static double Clamp(double x, double min, double max)
        {
            if (x < min)
                x = min;
            else if (x > max)
                x = max;

            return x;
        }
            
    }
}