using System;
using MGSharp.Core.MGCellModels;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGModels;
using System.Threading.Tasks;
using System.Collections.Generic;
using MGSharp.Core.Helpers;
using System.Text;
using System.Linq;

namespace MGSharp.Core.MGCellPopulation
{
    class CellPopulation
    {
        public int populationSize;
        public int maxPopulationSize;
        public int lastFramePopulationSize;
        public List<Tissue> tissues;
        public MGCell[] cells;
        public int[] sigma;
        public Vector dim;
        public Vector halfDim;
        int appliedForces;
        int frame = 0;
        StringBuilder interCellEdges;
        public Vector reference = new Vector();
        List<Vector> subset;
        public string name;
        List<MGCell> mGCells;
        public int[] contacts;
        public int[] contactsIndices;

        public CellPopulation(){}

        public CellPopulation(int popSize, int popMaxSize)
        {
            //Console.WriteLine("Creating Population ...");
            CreatePopulation(popSize, popMaxSize);
        }

        public CellPopulation(int popSize, Vector dim)
        {
            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = (int)(dim.x * dim.y * dim.z);
            cells = new MGCell[maxPopulationSize];
            sigma = new int[maxPopulationSize];
            this.dim = dim;
            halfDim = dim / 2;

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];
            for (int i = 0; i < populationSize; i++)
            {
                cells[i] = new MGCell(i);
                contactsIndices[i] = i;
                contacts[i] = 0;
            }
        }

        public CellPopulation(int popSize, int popMaxSize, List<Tissue> tissueList)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];

            tissues = tissueList;
            int start = 0;
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    //int k = i == 0 ? 0 : i - 1;
                    //int index = i * tissueList[k].populationSize + j;

                    int index = start + j;
                    cells[index] = new MGCell(index, tissueList[i].mesh);
                    cells[index].tissueName = tissueList[i].name;
                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;

                    contactsIndices[index] = index;
                    contacts[index] = 0;
                }
            }
        }

        public CellPopulation(int popSize, int popMaxSize, List<Tissue> tissueList, short cellType)
        {
            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];

            tissues = tissueList;
            int start = 0;
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    int index = start + j;
                    switch (cellType)
                    {
                        case 1: cells[index] = new RedBloodCell(index, tissueList[i].mesh); break;
                        case 2: cells[index] = new EpithelialCell(index, tissueList[i].mesh); break;
                        default: cells[index] = new MGCell(index, tissueList[i].mesh); break;
                    }
                    //cells[index] = new MGCell(index, tissueList[i].mesh);
                    cells[index].tissueName = tissueList[i].name;
                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;

                    contactsIndices[index] = index;
                    contacts[index] = 0;
                }
            }
        }

        public CellPopulation(int popSize, int popMaxSize, List<Tissue> tissueList, short cellType, bool randomCycleTime)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            //Console.WriteLine(popMaxSize);
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];

            tissues = tissueList;
            int start = 0;
            Random rand = new Random(0);
            
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    int index = start + j;

                    cells[index] = new MGCell(index, tissueList[i].mesh);
                    cells[index].tissueId = i;
                    cells[index].tissueName = tissueList[i].name;

                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                    cells[index].cycleTime = 0;// rand.Next(MGModel.cellCyclePeriod);

                    contactsIndices[index] = index;
                    contacts[index] = 0;
                }
            }
        }

        public CellPopulation(int popMaxSize, List<Tissue> tissueList)
        {
            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = 0;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[popMaxSize];
            double d = Math.Round(Math.Pow(maxPopulationSize, 1f / 3f));
            dim = new Vector(d, d, d);
            halfDim = dim / 2;

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];

            tissues = tissueList;
            int start = 0;
            for (int i = 0; i < tissueList.Count; i++)
            {
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissueList[i - 1].populationSize;
                for (int j = 0; j < tissueList[i].populationSize; j++)
                {
                    int index = start + j;
                    cells[index] = new MGCell(index, i, tissueList[i].cells[j]);
                    cells[index].tissueName = tissueList[i].name;
                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;

                    contactsIndices[index] = index;
                    contacts[index] = 0;
                }
                populationSize += tissueList[i].populationSize;
            }
        }

        public CellPopulation(string folder)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            CreateCellPopulation(folder);
        }
        
        public void CreatePopulation(int popSize, int popMaxSize)
        {
            //interCellEdges = new StringBuilder();
            //interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            Simulator.interCellEdges = new StringBuilder();
            Simulator.interCellEdges.AppendLine("frame,cell1,p1,cell2,p2");

            populationSize = popSize;
            maxPopulationSize = popMaxSize;
            cells = new MGCell[popMaxSize];
            sigma = new int[maxPopulationSize];
            dim = new Vector(Math.Pow(maxPopulationSize, 1f / 3f), Math.Pow(maxPopulationSize, 1f / 3f), Math.Pow(maxPopulationSize, 1f / 3f));
            halfDim = dim / 2;

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];

            for (int i = 0; i < populationSize; i++)
            {
                cells[i] = new MGCell(i);
                //Console.WriteLine("Created cell " + i + ", " + cells[i].cellId);
            }

            tissues = new List<Tissue>();
        }
        
        void CreateCellPopulation(string folder)
        {
            CellProvider provider = new CellProvider();
            provider.ReadParameters(folder);
            provider.ReadResults();
            provider.ReadParticleNeighbourhoods();
            provider.PlayResults(0);

            populationSize = provider.popSize;
            maxPopulationSize = provider.maxPopSize;

            cells = new MGCell[provider.maxPopSize];
            sigma = new int[provider.maxPopSize];

            contacts = new int[maxPopulationSize];
            contactsIndices = new int[maxPopulationSize];

            int start = 0;
            
            tissues = new List<Tissue>(provider.nbOfCellTypes);
            
            Random rand = new Random(0);
            for (int i = 0; i < provider.nbOfCellTypesPerFrame[provider.nbOfFrames-1]; i++)
            {
                //Console.WriteLine(provider.nbOfCellsPerTypePerFrame[provider.nbOfFrames - 1][i]);
                tissues.Add(new Tissue(provider.nbOfCellsPerTypePerFrame[provider.nbOfFrames - 1][i], provider.maxPopSize));
                tissues[i].indices = new int[tissues[i].maxPopulationSize];

                start += i == 0 ? 0 : tissues[i - 1].populationSize;
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    int index = start + j;

                    //Console.WriteLine(index);
                    cells[index] = new MGCell(index, i, provider.meshPerCell[index]);
                    //Console.WriteLine(index + ", " + provider.meshPerCell[index].vertexCount());
                    cells[index].tissueId = i;

                    tissues[i].indices[j] = index;
                    tissues[i].cells[j] = cells[index];
                    sigma[index] = index;
                    cells[index].cycleTime = rand.Next(MGModel.cellCyclePeriod);
                    tissues[i].reference += cells[index].ComputeCentreFromMesh();
                    //Console.WriteLine(cells[index].cycleTime);

                    contactsIndices[index] = index;
                    contacts[index] = 0;

                }
                tissues[i].mesh = provider.meshPerCell[start];
                tissues[i].reference /= tissues[i].populationSize;
            }

            for(int i=0; i<provider.externalEdges.Count; i++)
            {
                Edge edge = new Edge(cells[(int)provider.externalEdges[i].x].vertices[(int)provider.externalEdges[i].y],
                                        cells[(int)provider.externalEdges[i].z].vertices[(int)provider.externalEdges[i].w]);

                cells[(int)provider.externalEdges[i].x].externalEdges.add(edge);
                cells[(int)provider.externalEdges[i].z].externalEdges.add(new Edge(edge.ends[1], edge.ends[0]));

                cells[(int)provider.externalEdges[i].x].vertices[(int)provider.externalEdges[i].y].externalNeighbours.Add(new int[] {
                        (int)provider.externalEdges[i].z,
                        (int)provider.externalEdges[i].w
                });

                cells[(int)provider.externalEdges[i].z].vertices[(int)provider.externalEdges[i].w].externalNeighbours.Add(new int[] {
                        (int)provider.externalEdges[i].x,
                        (int)provider.externalEdges[i].y
                });

                if((int)provider.externalEdges[i].x == 0 && (int)provider.externalEdges[i].y == 26)
                {
                    //Console.WriteLine(provider.externalEdges[i].z + ", " + provider.externalEdges[i].w);
                }
            }

            for(int i=0; i<populationSize; i++)
            {
                contacts[i] = cells[i].externalEdges.getCount();
            }
        }

        public void AddCell(int tissueId, MGCell cell)
        {
            if (populationSize < maxPopulationSize)
            {
                tissues[tissueId].AddCell(cell);
                cells[populationSize] = cell;
                populationSize++;
            }
        }

        public void AddTissue(Tissue tissue)
        {
            int i = tissues.Count;
            tissue.indices = new int[tissue.maxPopulationSize];

            int start = populationSize;
            for (int j = 0; j < tissue.populationSize; j++)
            {
                int index = start + j;
                cells[index] = new MGCell(index, tissue.mesh);
                
                tissue.indices[j] = index;
                tissue.cells[j] = cells[index];
                sigma[index] = index;
            }
            tissues.Add(tissue);
            populationSize += tissue.populationSize;
        }

        public void PositionCells(bool random)
        {
            Vector centre = new Vector(0,0,0);

            if (!random)
            {
                for (int i = 0; i < populationSize; ++i)
                {
                    double x = Helper.NCubedInN[i].x * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    double y = Helper.NCubedInN[i].y * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    double z = Helper.NCubedInN[i].z * MGModel.Rcell * cells[i].innerRadius.z * 2;

                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(centre);

                    reference.y = y;
                }
            }
            else
            {
                for (int i = 0; i < populationSize; ++i)
                {
                    Random rnd = new Random();
                    double x = rnd.Next(0, (int)dim.x) * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    double y = rnd.Next(0, (int)dim.y) * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    double z = rnd.Next(0, (int)dim.z) * MGModel.Rcell * cells[i].innerRadius.z * 2;
                    
                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(centre);
                }
            }
        }

        public void PositionCells(bool random, Vector refer)
        {
            Vector centre = new Vector(0, 0, 0);
            reference = refer;

            if (!random)
            {
                double x=0, y=0, z = 0;
                for (int i = 0; i < populationSize; ++i)
                {
                    x = Helper.NCubedInN[i].x * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    y = Helper.NCubedInN[i].y * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    z = Helper.NCubedInN[i].z * MGModel.Rcell * cells[i].innerRadius.z * 2;

                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(centre+reference);
                }
                reference.y += y;
            }
            else
            {
                for (int i = 0; i < populationSize; ++i)
                {
                    Random rnd = new Random();
                    double x = rnd.Next(0, (int)dim.x) * MGModel.Rcell * cells[i].innerRadius.x * 2;
                    double y = rnd.Next(0, (int)dim.y) * MGModel.Rcell * cells[i].innerRadius.y * 2;
                    double z = rnd.Next(0, (int)dim.z) * MGModel.Rcell * cells[i].innerRadius.z * 2;

                    centre.x = x - halfDim.x;
                    centre.y = y - halfDim.y;
                    centre.z = z - halfDim.z;

                    cells[i].PositionCell(reference + centre);
                }
            }
        }
        
        public void SetCellNeighbourhood(int cellId)
        {
            MGCell cell = cells[cellId];
            cell.neighbours.Clear();

            if (MGModel.hexagonalNeighbourhoodForCells)
            {
                cell.neighbours = Helper.GetHexagonalNeighbourhood3D(cellId, maxPopulationSize, populationSize);
                //Console.WriteLine("Hexagonal");
            }
            else if (MGModel.mooreNeighbourhoodForCells)
            {
                //cell.neighbours = Helper.GetMooreNeighbourhood3D(cellId, maxPopulationSize, populationSize);
                cell.neighbours = Helper.GetMooreNeighbourhood3D(cellId, dim, populationSize);
                //Console.WriteLine(maxPopulationSize + ", " + populationSize);
            }
            else
            {
                //Console.WriteLine("Population Size: " + populationSize);
                //Console.WriteLine("Cell ID: " + cellId);
                for (int i = 0; i < populationSize; i++)
                {
                    double distance = Vector.Distance(cell.GetPosition(), cells[i].GetPosition());

                    if (i != cellId && distance < MGModel.maximumNeighbourDistance && cell.neighbours.Count < 26)
                    {
                        cell.neighbours.Add(i);
                        //Console.WriteLine(cellId + ":" + cell.GetPosition() + ", " + i + ":" + cells[i].GetPosition() + ", " + distance + ", " + MGModel.maximumNeighbourDistance);
                    }
                }
            }
        }
        
        public void DynamiseSystem()
        {
            Parallel.For(0, populationSize, i =>
            {
                cells[i].ComputeCentreFromMesh();
                for(int j=0; j < cells[i].vertexCount(); j++)
                {
                    cells[i].vertices[j].externalForces = new Vector();
                }
            });

            if (MGModel.staticNeighbourhood)
            {
                if (Simulator.frame < 1 || Simulator.frame == MGModel.nextNeighboursSearchFrame)
                {
                    if (!MGModel.readNeigbours)
                    {
                        for (int i = 0; i < populationSize; i++)
                        {
                            SetCellNeighbourhood(i);
                            cells[i].ComputeForces();
                            SetInterCellEdges(i);
                            //SetInterCellEdgesWithTriangles(i);
                        }
                    }
                    else
                    {
                        //Console.WriteLine("I read Neighbourhoods !");
                        for (int i = 0; i < populationSize; i++)
                        {
                            //SetCellNeighbourhood(i);
                            //Console.WriteLine(i + ", " + cells[i].externalEdges.getCount());
                            cells[i].ComputeForces();
                        }
                    }
                }

                for (int i = 0; i < populationSize; i++)
                {
                    cells[i].ComputeForces();
                }

                for (int i = 0; i < populationSize; i++)
                {
                    cells[i].ComputeExternalForces();
                }

                if (Simulator.frame == 222000)
                {
                    CloseEPI2();
                    for (int i = 0; i < populationSize; i++)
                    {
                        cells[i].ComputeExternalForces();
                    }
                }

                for (int i = 0; i < populationSize; i++)
                {
                    cells[i].Dynamise();
                }
            }
            else
            {
                MGModel.nextNeighboursSearchFrame = Simulator.frame;

                Parallel.For(0, populationSize, i =>
                {
                    cells[i].ComputeForces();
                    if (Simulator.frame == MGModel.nextNeighboursSearchFrame)
                    {
                        SetCellNeighbourhood(i);
                        SetInterCellEdges(i);
                        //SetInterCellEdgesWithTriangles(i);
                    }
                });
                
                
                Parallel.For(0, populationSize, i =>
                {
                    cells[i].ComputeExternalForces();
                });

                Parallel.For(0, populationSize, i =>
                {
                    cells[i].Dynamise();
                });
            }
            frame++;
        }

        public void SetInterCellEdges(int cellId)
        {
            Edge edge;
            double distance = 0;
            cells[cellId].externalEdges.clear();

            // Permutation of particles
            sigma = new int[cells[cellId].nbOfParticles];
            
            if (Simulator.frame % 10 == 0)
            {
                Randomize();
            }

            for (int i = 0; i < cells[cellId].nbOfParticles; i++)
            {
                double minDistance = 1000;

                // To make sure all external neighbours have their external forces defined
                cells[cellId].vertices[i].externalForces = new Vector();
                cells[cellId].vertices[i].externalNeighbours = new List<int[]>();

                int countN = 0;
                for (int j = 0; j < cells[cellId].neighbours.Count; j++)
                {
                    int l = cells[cellId].neighbours[j];
                    for (int k = 0; k < cells[l].nbOfParticles; k++)
                    {
                        distance = Vector.Distance(cells[cellId].vertices[i].GetPosition(), cells[l].vertices[k].GetPosition());

                        if (distance <= minDistance)
                        {
                            minDistance = distance;
                            cells[cellId].vertices[i].closestNeighbour = cells[l].vertices[k].Clone();
                            //cells[cellId].vertices[i].closestNeighbour.v = new Vector(cells[l].vertices[k].v);
                        }
                          

                        if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId] && countN < MGModel.maxExternalNeighbours)
                        //distance = Vector.Distance(cells[cellId].vertices[i].GetPosition(), cells[l].vertices[sigma[k]].GetPosition());
                        //if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId] && countN < MGModel.maxExternalNeighbours)
                        {
                            countN++;
                            //Console.WriteLine(distance);
                            edge = new Edge(cells[cellId].vertices[i], cells[l].vertices[k]);
                            //edge = new Edge(cells[cellId].vertices[i], cells[l].vertices[sigma[k]]);
                            cells[cellId].externalEdges.add(edge);
                            cells[cellId].vertices[i].externalNeighbours.Add(new int[] { l, k });

                            cells[l].externalEdges.add(new Edge(edge.ends[1], edge.ends[0]));
                            cells[l].vertices[k].externalNeighbours.Add(new int[] { cellId, i });
                        }
                    }
                }
            }

            contacts[cellId] = cells[cellId].neighbours.Count;
            contactsIndices[cellId] = cellId;
        }

        public void SetInterCellEdgesWithTriangles(int cellId)
        {
            //Console.WriteLine("Searching ...");
            cells[cellId].externalEdges.clear();

            double distance = 0;

            // Permutation of particles
            //sigma = new int[cells[cellId].nbOfParticles];

            for (int i = 0; i < cells[cellId].nbOfParticles; i++)
            {
                // To make sure all external neighbours have their external forces defined
                cells[cellId].vertices[i].externalForces = new Vector();
                cells[cellId].vertices[i].externalNeighbours = new List<int[]>();

                for (int j = 0; j < cells[cellId].neighbours.Count; j++)
                {
                    //int closestTriangle;
                    int l = cells[cellId].neighbours[j];

                    for (int k = 0; k < cells[l].faceCount(); k++)
                    {
                        distance = cells[l].faces[k].DistanceToTriangle(cells[cellId].vertices[i].v);
                        if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId])
                        {
                            //Console.WriteLine(cellId + ", " + i + ", " + l + ", " + k + ", " + distance);
                            Vector P = new Vector();

                            int key = 0;
                            bool test = Mesh.RayIntersectsTriangle(key, cells[cellId].vertices[i].v,
                                                    cells[cellId].vertices[i].force * 100,
                                                    cells[l].faces[k],
                                                    out P);
                            if (!test)
                            {
                                test = Mesh.RayIntersectsTriangle(key, cells[cellId].vertices[i].v,
                                                    -cells[cellId].vertices[i].force * 100,
                                                    cells[l].faces[k],
                                                    out P);
                            }
                                                    
                            //Console.WriteLine(P);
                            //distance = Vector.Distance(cells[cellId].vertices[i].v, P);

                            if (test)
                            {
                                Vertex A = cells[l].faces[k].vertices[0];
                                Vertex B = cells[l].faces[k].vertices[1];
                                Vertex C = cells[l].faces[k].vertices[2];

                                //double Aabc = ((A.v - B.v) ^ (A.v - C.v)).norm() / 2;
                                double gammaA = ((P - B.v) ^ (P - C.v)).norm() / 2;
                                double gammaB = ((P - A.v) ^ (P - C.v)).norm() / 2;
                                double gammaC = ((P - A.v) ^ (P - B.v)).norm() / 2;

                                double maxGamma = 0;
                                if (gammaA >= gammaB && gammaA >= gammaC)
                                    maxGamma = gammaA;
                                if (gammaB >= gammaA && gammaB >= gammaC)
                                    maxGamma = gammaB;
                                if (gammaC >= gammaA && gammaC >= gammaB)
                                    maxGamma = gammaC;

                                gammaA /= maxGamma;
                                gammaB /= maxGamma;
                                gammaC /= maxGamma;

                                //Console.WriteLine("Tests passed !");
                                //Console.WriteLine(gammaA + ", " + gammaB + ", " + gammaC);

                                Edge edgeA = new Edge(cells[cellId].vertices[i], A, gammaA);
                                Edge edgeB = new Edge(cells[cellId].vertices[i], B, gammaB);
                                Edge edgeC = new Edge(cells[cellId].vertices[i], C, gammaC);

                                cells[cellId].externalEdges.add(edgeA);
                                cells[cellId].externalEdges.add(edgeB);
                                cells[cellId].externalEdges.add(edgeC);

                                cells[cellId].vertices[i].externalNeighbours.Add(new int[] { A.cellId, A.pos });
                                cells[cellId].vertices[i].externalNeighbours.Add(new int[] { B.cellId, B.pos });
                                cells[cellId].vertices[i].externalNeighbours.Add(new int[] { C.cellId, C.pos });

                                //n = cells[l1].vertices[k1].faces.getCount();
                            }
                            //n++;

                        }
                    }
                }
            }
        }

        public void SetInterCellEdgesIfNotExisting(int cellId)
        {
            Edge edge;
            double distance = 0;
            //cells[cellId].externalEdges.clear();

            for (int i = 0; i < cells[cellId].nbOfParticles; i++)
            {
                // To make sure all external neighbours have their external forces defined
                cells[cellId].vertices[i].externalForces = new Vector();

                if (cells[cellId].vertices[i].externalNeighbours.Count == 0)
                {
                    int countN = 0;
                    for (int j = 0; j < cells[cellId].neighbours.Count; j++)
                    {
                        int l = cells[cellId].neighbours[j];
                        for (int k = 0; k < cells[l].nbOfParticles; k++)
                        {
                            distance = Vector.Distance(cells[cellId].vertices[i].GetPosition(), cells[l].vertices[k].GetPosition());
                            if (distance < MGModel.DInteraction[cells[cellId].tissueId, cells[l].tissueId] && countN < MGModel.maxExternalNeighbours)
                            {
                                countN++;

                                edge = new Edge(cells[cellId].vertices[i], cells[l].vertices[k]);
                                cells[cellId].externalEdges.add(edge);
                                cells[cellId].vertices[i].externalNeighbours.Add(new int[] { l, k });
                                //Console.WriteLine("Woo hoo");

                                //interCellEdges.AppendLine( frame + "," + cellId + "," + i + "," + l + "," + k);
                            }
                        }
                    }
                }
            }
        }

        public void CloseEPI2()
        {
            Cylinder34 cy = new Cylinder34(new Vector(.5f, 1f, .5f), false);
            int[] lookUp = cy.LookUp();

            int[] ring = new int[] { 58, 59, 60, 61, 62, 65, 66, 67, 68, 69, 72, 73, 74, 75, 76, 79, 80, 81, 82, 83, 86, 87, 88, 89, 90 };
            /*
            List<int> TE = new List<int>(){ 50, 51, 52, 53, 54,
                                            55, 56, 57, 58, 59,
                                            60, 61, 62, 63, 64,
                                            65, 66, 67, 68, 69,
                                            70, 71, 72, 73, 74};
            //*/

            //*
            List<int> TE = new List<int>();
            for (int i = 0; i < 49; i++)
            {
                TE.Add(50 + i);
            }
            //*/

            int[] pFoyer = { 6, 7, 8, 11, 12, 13, 16, 17, 18 };
            int[] pZmin = { 0, 1, 2, 3, 4};
            int[] pXmax = { 4, 9, 14, 19, 24 };
            int[] pZmax = { 24, 23, 22, 21, 20 };
            int[] pXmin = { 20, 15, 10, 5, 0 };

            //int[] cSup = { 25, 27, 29, 31, 32 };
            int[] cSup = { 24, 25, 26, 27, 28, 29, 30, 31, 32 };
            int[] cZmin = cy.CZmin().ToArray();
            int[] cZmax = cy.CZmax().ToArray();
            int[] cXmin = cy.CXmin().ToArray();
            int[] cXmax = cy.CXmax().ToArray();

            int edgeCount = 0;
            List<Edge> newEdges = new List<Edge>();
            for (int i = 0; i < pZmin.Length; i++)
            {
                for (int j = 0; j < cZmin.Length; j++)
                {
                    //Console.WriteLine(pZmin[i] + ", " + cZmin[j]);
                    Edge edge = new Edge(cells[pZmin[i]].vertices[cZmin[j]], cells[pZmin[i] + 25].vertices[lookUp[cZmin[j]]]);
                    newEdges.Add(edge);
                    edgeCount++;
                }   
            }

            for (int i = 0; i < pZmax.Length; i++)
            {
                for (int j = 0; j < cZmax.Length; j++)
                {
                    Edge edge = new Edge(cells[pZmax[i]].vertices[cZmax[j]], cells[pZmax[i] + 25].vertices[lookUp[cZmax[j]]]);
                    //Edge edge = new Edge(cells[pZmax[i]].vertices[cZmax[j]], cells[pZmax[i] + 25].vertices[cZmax[j]]);
                    newEdges.Add(edge);
                    edgeCount++;
                }
            }

            for (int i = 0; i < pXmin.Length; i++)
            {
                for (int j = 0; j < cXmin.Length; j++)
                {
                    Edge edge = new Edge(cells[pXmin[i]].vertices[cXmin[j]], cells[pXmin[i] + 25].vertices[lookUp[cXmin[j]]]);
                    //Edge edge = new Edge(cells[pZmin[i]].vertices[cXmin[j]], cells[pZmin[i] + 25].vertices[cXmin[j]]);
                    newEdges.Add(edge);
                    edgeCount++;
                }
            }

            for (int i = 0; i < pXmax.Length; i++)
            {
                for (int j = 0; j < cXmax.Length; j++)
                {
                    Edge edge = new Edge(cells[pXmax[i]].vertices[cXmax[j]], cells[pXmax[i] + 25].vertices[lookUp[cXmax[j]]]);
                    newEdges.Add(edge);
                    edgeCount++;
                }
            }

            for (int i = 0; i < pFoyer.Length; i++)
            {
                for (int j = 0; j < cSup.Length; j++)
                {
                    Edge edge = new Edge(cells[pFoyer[i]].vertices[cSup[j]], cells[pFoyer[i] + 25].vertices[lookUp[cSup[j]]]);

                    for(int n=0; n<edge.ends[1].externalNeighbours.Count; n++)
                    {
                        int l = edge.ends[1].externalNeighbours[n][0];
                        int k = edge.ends[1].externalNeighbours[n][1];

                        Edge e1 = new Edge(edge.ends[0], cells[l].vertices[k]);
                        newEdges.Add(e1);
                    }

                    newEdges.Add(edge);
                    edgeCount++;
                }
            }

            //newEdges = newEdges.Distinct().ToList<Edge>();

            newEdges = Helper.Distinct(newEdges);
            int count = 0;
            //*
            for (int i = 0; i < newEdges.Count; i++)
            {
                count++;
                Vector midPoint = (newEdges[i].ends[0].v + newEdges[i].ends[1].v) * .5;

                newEdges[i].ends[0].v = new Vector(midPoint);
                newEdges[i].ends[1].v = new Vector(midPoint);

                Vertex v0 = newEdges[i].ends[0];
                Vertex v1 = newEdges[i].ends[1];

                /*
                if(v0.cellId == 22 && v0.pos == 4 && v1.cellId == 47 && v1.pos == 2)
                {
                    for (int j = 0; j < v1.externalNeighbours.Count; j++)
                    {
                        int l = v1.externalNeighbours[j][0];
                        int k = v1.externalNeighbours[j][1];

                        //Console.WriteLine(l + ", " + k);
                    }
                }
                */    

                //*
                for (int j = 0; j < v0.externalNeighbours.Count; j++)
                {
                    int l = v0.externalNeighbours[j][0];
                    int k = v0.externalNeighbours[j][1];

                    if(l != v1.cellId)
                    {
                        Edge newEdge = new Edge(v1, cells[l].vertices[k]);
                        cells[v1.cellId].externalEdges.add(newEdge);
                        cells[l].externalEdges.add(new Edge(cells[l].vertices[k], v1));

                        //v1.externalNeighbours.Add(new int[] { l, k });
                    }
                }

                for (int j = 0; j < v1.externalNeighbours.Count; j++)
                {
                    int l = v1.externalNeighbours[j][0];
                    int k = v1.externalNeighbours[j][1];

                    //cells[l].vertices[k].v = new Vector(v1.v);

                    if (l != v0.cellId)
                    {
                        if (TE.Contains(l))// && !TE.Contains(l))
                        {
                            Edge newEdge = new Edge(v1, cells[l].vertices[k]);
                            cells[v1.cellId].externalEdges.remove(newEdge);
                            cells[l].externalEdges.remove(new Edge(cells[l].vertices[k], v1));
                        }
                        else
                        {
                            Edge newEdge = new Edge(v0, cells[l].vertices[k]);
                            cells[v0.cellId].externalEdges.add(newEdge);
                            cells[l].externalEdges.add(new Edge(cells[l].vertices[k], v0));

                            //v0.externalNeighbours.Add(new int[] { l, k });
                        }
                    }
                }
                //*/

                v0.externalNeighbours.Add(new int[] { v1.cellId, v1.pos });
                v1.externalNeighbours.Add(new int[] { v0.cellId, v0.pos });

                cells[newEdges[i].ends[0].cellId].externalEdges.add(newEdges[i]);
                cells[newEdges[i].ends[1].cellId].externalEdges.add(new Edge(newEdges[i].ends[1], newEdges[i].ends[0]));
            }
            //*/

            int[] corners = new int[] { 0, 4, 24, 20};
            for(int i=0; i<corners.Length; i++)
            {
                Cylinder34 cl = new Cylinder34(new Vector(.5f, 1f, .5f), false);
                cl.ApicalConstrictionWithPositionChange(.4f, Vector.up);

                Cylinder34 clDown = new Cylinder34(new Vector(.5f, 1f, .5f), false);
                clDown.ApicalConstrictionWithPositionChange(.4f, Vector.down);

                Vector normal = new Vector();
                Vector normal25 = new Vector();
                if (i == 0)
                {
                    normal = new Vector(cl.vertices[0].v.x, 0, cl.vertices[0].v.z);
                    //normal25 = new Vector(cl.vertices[4].v.x, 0, cl.vertices[4].v.z);
                    normal25 = new Vector(cl.vertices[6].v.x, 0, cl.vertices[6].v.z);
                }
                if (i == 1)
                {
                    normal = new Vector(cl.vertices[2].v.x, 0, cl.vertices[2].v.z);
                    //normal25 = new Vector(cl.vertices[6].v.x, 0, cl.vertices[6].v.z);
                    normal25 = new Vector(cl.vertices[4].v.x, 0, cl.vertices[4].v.z);
                }
                if (i == 2)
                {
                    normal = new Vector(cl.vertices[4].v.x, 0, cl.vertices[4].v.z);
                    //normal25 = new Vector(cl.vertices[0].v.x, 0, cl.vertices[0].v.z);
                    normal25 = new Vector(cl.vertices[2].v.x, 0, cl.vertices[2].v.z);
                } 
                if (i == 3)
                {
                    normal = new Vector(cl.vertices[6].v.x, 0, cl.vertices[6].v.z);
                    //normal25 = new Vector(cl.vertices[2].v.x, 0, cl.vertices[2].v.z);
                    normal25 = new Vector(cl.vertices[0].v.x, 0, cl.vertices[0].v.z);
                }

                cells[corners[i]].SetElengths(cl.Slice(normal, true));
                cells[corners[i] + 25].SetElengths(clDown.Slice(normal25, false));

                //for(int )

                //cells[corners[i]].SetElengths(cl.Slice(normal, true));
                //cells[corners[i] + 25].SetElengths(clDown.Slice(normal25, true));

                //cells[corners[i] + 25].SetEdgeELengths();
            }
        }

        public void StochasticDynamiseSystem()
        {
            Parallel.For(0, populationSize, i =>
            {
                cells[i].ResetCell();
                SetCellNeighbourhood(i);
                SetInterCellEdges(i);
                cells[i].ChangeShape();
                cells[i].ComputeForces();
                cells[i].Dynamise();
            });
        }

        public void Cleavage(bool staticShape)
        {
            Parallel.For(0, populationSize, i =>
            {
                cells[i].ComputeCentreFromMesh();
                for (int j = 0; j < cells[i].vertexCount(); j++)
                {
                    cells[i].vertices[j].externalForces = new Vector();
                    cells[i].vertices[j].globalForces = new Vector();
                }
            });

            Vector normal = populationSize == 1 ? new Vector(1, 0, 0) : Vector.zero;
            for (int i = 0; i < populationSize; i++)
            {
                cells[i].CellCycle(appliedForces, normal, staticShape);
            }
        }

        public void Randomize()
        {
            Randomizer.Randomize<int>(sigma);
        }

        public void Translate(Vector u)
        {
            for(int i = 0; i < populationSize; i++)
            {
                cells[i].Translate(u);
            }
        }

        public void LogNumbers()
        {
            PersistantNumbers pv = 
                new PersistantNumbers(Simulator.frame, tissues.Count, populationSize);
            Simulator.numbers.Add(pv);
            for(int i = 0; i < tissues.Count; i++)
            {
                pv = new PersistantNumbers(Simulator.frame, tissues[i].id, tissues[i].populationSize);
                Simulator.numbers.Add(pv);
                for (int j=0; j < tissues[i].populationSize; j++)
                {
                    pv = new PersistantNumbers(Simulator.frame, tissues[i].cells[j].cellId, tissues[i].cells[j].nbOfParticles, i);
                    Simulator.numbers.Add(pv);
                }
            }
        }

        public void LogFinalNumbers()
        {
            PersistantNumbers pv =
                new PersistantNumbers(Simulator.frame, tissues.Count, populationSize);
            Simulator.finalNumbers.Add(pv);
            for (int i = 0; i < tissues.Count; i++)
            {
                pv = new PersistantNumbers(Simulator.frame, tissues[i].id, tissues[i].populationSize);
                Simulator.finalNumbers.Add(pv);
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    pv = new PersistantNumbers(Simulator.frame, tissues[i].cells[j].cellId, tissues[i].cells[j].nbOfParticles, i);
                    Simulator.finalNumbers.Add(pv);
                }
            }
        }
        
        public void LogMetrics()
        {
            PersistantMetrics pm = new PersistantMetrics(Simulator.frame, Metrics.ElasticEnergy());//Metrics.MouseEmbryo()); //Metrics.ElasticEnergy());//Metrics.ForcesRatio());

            Simulator.metrics.Add(pm);
        }

        public void LogNeighbourhoods()
        {
            for(int i=0; i<populationSize; i++)
            {
                for(int j=0; j<cells[i].externalEdges.getCount(); j++)
                {
                    Simulator.interCellEdges.AppendLine(Simulator.frame + ";" + 
                                                cells[i].externalEdges[j].ends[0].cellId + ";" +
                                                cells[i].externalEdges[j].ends[0].pos + ";" +
                                                cells[i].externalEdges[j].ends[1].cellId + ";" +
                                                cells[i].externalEdges[j].ends[1].pos + ";" +
                                                cells[i].externalEdges[j].gamma);
                }
            }
        }

        public void LogState()
        {
            PersistantVertex pv;
            for (int i = 0; i < tissues.Count; i++)
            {
                for(int j=0; j<tissues[i].populationSize; j++)
                {
                    for(int k=0; k<tissues[i].cells[j].nbOfParticles; k++)
                    {
                        pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].vertices[k].id, tissues[i].cells[j].cellId, tissues[i].cells[j].vertices[k].v);
                        Simulator.states.Add(pv.Clone());
                    }
                    pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].nuclei[0].id, tissues[i].cells[j].cellId, tissues[i].cells[j].nuclei[0].v);
                    Simulator.states.Add(pv.Clone());
                }
            }
        }

        public void LogFinalState()
        {
            PersistantVertex pv;
            for (int i = 0; i < tissues.Count; i++)
            {
                for (int j = 0; j < tissues[i].populationSize; j++)
                {
                    for (int k = 0; k < tissues[i].cells[j].nbOfParticles; k++)
                    {
                        pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].vertices[k].id, tissues[i].cells[j].cellId, tissues[i].cells[j].vertices[k].v);
                        Simulator.finalStates.Add(pv.Clone());
                    }
                    pv = new PersistantVertex(Simulator.frame, tissues[i].cells[j].nuclei[0].id, tissues[i].cells[j].cellId, tissues[i].cells[j].nuclei[0].v);
                    Simulator.finalStates.Add(pv.Clone());
                }
            }
        }
        //*/
    }
}