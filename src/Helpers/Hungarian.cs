﻿using System;

namespace MGSharp.Core.Helpers
{
    /// <summary>
    /// Hungarian Algorithm.
    /// </summary>
    public static class HungarianAlgorithm
    {
        /// <summary>
        /// Finds the optimal assignments for a given matrix of agents and costed tasks such that the total cost is minimized.
        /// </summary>
        /// <param name="costs">A cost matrix; the element at row <em>i</em> and column <em>j</em> represents the cost of agent <em>i</em> performing task <em>j</em>.</param>
        /// <returns>A matrix of assignments; the value of element <em>i</em> is the column of the task assigned to agent <em>i</em>.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="costs"/> is null.</exception>
        public static int[] FindAssignments(this float[,] costs)
        {
            if (costs == null)
                throw new ArgumentNullException("Empty");

            var h = costs.GetLength(0);
            var w = costs.GetLength(1);

            for (var i = 0; i < h; i++)
            {
                float min = int.MaxValue;

                for (var j = 0; j < w; j++)
                {
                    min = Math.Min(min, costs[i, j]);
                }

                for (var j = 0; j < w; j++)
                {
                    costs[i, j] -= min;
                }
            }

            var masks = new byte[h, w];
            var rowsCovered = new bool[h];
            var colsCovered = new bool[w];

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (costs[i, j] == 0 && !rowsCovered[i] && !colsCovered[j])
                    {
                        masks[i, j] = 1;
                        rowsCovered[i] = true;
                        colsCovered[j] = true;
                    }
                }
            }

            HungarianAlgorithm.ClearCovers(rowsCovered, colsCovered, w, h);

            var path = new Location[w * h];
            var pathStart = default(Location);
            var step = 1;

            while (step != -1)
            {
                switch (step)
                {
                    case 1:
                        step = HungarianAlgorithm.RunStep1(masks, colsCovered, w, h);
                        break;

                    case 2:
                        step = HungarianAlgorithm.RunStep2(costs, masks, rowsCovered, colsCovered, w, h, ref pathStart);
                        break;

                    case 3:
                        step = HungarianAlgorithm.RunStep3(masks, rowsCovered, colsCovered, w, h, path, pathStart);
                        break;

                    case 4:
                        step = HungarianAlgorithm.RunStep4(costs, rowsCovered, colsCovered, w, h);
                        break;
                }
            }

            var agentsTasks = new int[h];

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (masks[i, j] == 1)
                    {
                        agentsTasks[i] = j;
                        break;
                    }
                }
            }

            return agentsTasks;
        }

        private static int RunStep1(byte[,] masks, bool[] colsCovered, int w, int h)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");
                //throw new ArgumentNullException(nameof(masks));

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (masks[i, j] == 1)
                        colsCovered[j] = true;
                }
            }

            var colsCoveredCount = 0;

            for (var j = 0; j < w; j++)
            {
                if (colsCovered[j])
                    colsCoveredCount++;
            }

            if (colsCoveredCount == h)
                return -1;

            return 2;
        }
        private static int RunStep2(float[,] costs, byte[,] masks, bool[] rowsCovered, bool[] colsCovered, int w, int h, ref Location pathStart)
        {
            if (costs == null)
                throw new ArgumentNullException("Empty");

            if (masks == null)
                throw new ArgumentNullException("Empty");

            if (rowsCovered == null)
                throw new ArgumentNullException("Empty");

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            while (true)
            {
                var loc = HungarianAlgorithm.FindZero(costs, rowsCovered, colsCovered, w, h);
                if (loc.row == -1)
                    return 4;

                masks[loc.row, loc.column] = 2;

                var starCol = HungarianAlgorithm.FindStarInRow(masks, w, loc.row);
                if (starCol != -1)
                {
                    rowsCovered[loc.row] = true;
                    colsCovered[starCol] = false;
                }
                else
                {
                    pathStart = loc;
                    return 3;
                }
            }
        }
        private static int RunStep3(byte[,] masks, bool[] rowsCovered, bool[] colsCovered, int w, int h, Location[] path, Location pathStart)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");

            if (rowsCovered == null)
                throw new ArgumentNullException("Empty");

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            var pathIndex = 0;
            path[0] = pathStart;

            while (true)
            {
                var row = HungarianAlgorithm.FindStarInColumn(masks, h, path[pathIndex].column);
                if (row == -1)
                    break;

                pathIndex++;
                path[pathIndex] = new Location(row, path[pathIndex - 1].column);

                var col = HungarianAlgorithm.FindPrimeInRow(masks, w, path[pathIndex].row);

                pathIndex++;
                path[pathIndex] = new Location(path[pathIndex - 1].row, col);
            }

            HungarianAlgorithm.ConvertPath(masks, path, pathIndex + 1);
            HungarianAlgorithm.ClearCovers(rowsCovered, colsCovered, w, h);
            HungarianAlgorithm.ClearPrimes(masks, w, h);

            return 1;
        }
        private static int RunStep4(float[,] costs, bool[] rowsCovered, bool[] colsCovered, int w, int h)
        {
            if (costs == null)
                throw new ArgumentNullException("Empty");

            if (rowsCovered == null)
                throw new ArgumentNullException("Empty");

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            var minValue = HungarianAlgorithm.FindMinimum(costs, rowsCovered, colsCovered, w, h);

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (rowsCovered[i])
                        costs[i, j] += minValue;
                    if (!colsCovered[j])
                        costs[i, j] -= minValue;
                }
            }
            return 2;
        }

        private static float FindMinimum(float[,] costs, bool[] rowsCovered, bool[] colsCovered, int w, int h)
        {
            if (costs == null)
                throw new ArgumentNullException("Empty");

            if (rowsCovered == null)
                throw new ArgumentNullException("Empty");

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            float minValue = int.MaxValue;

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (!rowsCovered[i] && !colsCovered[j])
                        minValue = Math.Min(minValue, costs[i, j]);
                }
            }

            return minValue;
        }
        private static int FindStarInRow(byte[,] masks, int w, int row)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");

            for (var j = 0; j < w; j++)
            {
                if (masks[row, j] == 1)
                    return j;
            }

            return -1;
        }
        private static int FindStarInColumn(byte[,] masks, int h, int col)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");

            for (var i = 0; i < h; i++)
            {
                if (masks[i, col] == 1)
                    return i;
            }

            return -1;
        }
        private static int FindPrimeInRow(byte[,] masks, int w, int row)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");

            for (var j = 0; j < w; j++)
            {
                if (masks[row, j] == 2)
                    return j;
            }

            return -1;
        }
        private static Location FindZero(float[,] costs, bool[] rowsCovered, bool[] colsCovered, int w, int h)
        {
            if (costs == null)
                throw new ArgumentNullException("Empty");

            if (rowsCovered == null)
                throw new ArgumentNullException("Empty");

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (costs[i, j] == 0 && !rowsCovered[i] && !colsCovered[j])
                        return new Location(i, j);
                }
            }

            return new Location(-1, -1);
        }
        private static void ConvertPath(byte[,] masks, Location[] path, int pathLength)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");

            if (path == null)
                throw new ArgumentNullException("Empty");

            for (var i = 0; i < pathLength; i++)
            {
                if (masks[path[i].row, path[i].column] == 1)
                {
                    masks[path[i].row, path[i].column] = 0;
                }
                else if (masks[path[i].row, path[i].column] == 2)
                {
                    masks[path[i].row, path[i].column] = 1;
                }
            }
        }
        private static void ClearPrimes(byte[,] masks, int w, int h)
        {
            if (masks == null)
                throw new ArgumentNullException("Empty");

            for (var i = 0; i < h; i++)
            {
                for (var j = 0; j < w; j++)
                {
                    if (masks[i, j] == 2)
                        masks[i, j] = 0;
                }
            }
        }
        private static void ClearCovers(bool[] rowsCovered, bool[] colsCovered, int w, int h)
        {
            if (rowsCovered == null)
                throw new ArgumentNullException("Empty");

            if (colsCovered == null)
                throw new ArgumentNullException("Empty");

            for (var i = 0; i < h; i++)
            {
                rowsCovered[i] = false;
            }

            for (var j = 0; j < w; j++)
            {
                colsCovered[j] = false;
            }
        }

        private struct Location
        {
            public readonly int row;
            public readonly int column;

            public Location(int row, int col)
            {
                this.row = row;
                this.column = col;
            }
        }
    }
}