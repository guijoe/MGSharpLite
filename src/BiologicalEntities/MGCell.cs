using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.Helpers;
using MGSharp.Core.LinearAlgebra;
using MGSharp.Core.MGModels;
//using Microsoft.Research.Oslo;

namespace MGSharp.Core.MGCellModels
{
    public class MGCell : Mesh
    {
        public int cellId;
        public int parent;
        public int tissueId;
        public int nbOfParticles;
        public Vertex[] nuclei;
        //public Vector centre;
        
        public List<int> neighbours;
        public EdgeSet externalEdges;
        public EdgeSet nucleusEdges;

        public Vector[] targetVertices;
        private Vector[] potentialNucleiPos;
        public float[,] particlesDistances;
        public int[] elongationAxis;

        public bool inGrowMode = true;
        public bool inDivisionMode = false;
        public bool hasDivided = false;
        public bool isNewBorn = true;

        public int cyclePeriod = 200;
        public int cycleTime = 0;

        public int appliedForces;
        public float Rcell;

        public float[] spins;
        public int[] sigma;
        public Vector polarisation;
        public Vector[] axis;

        public string tissueName;
        //public int tissueId;
        public int zeroNeighbours;
        public double divisionRate = 1;

        public double c;
        public double theta;
        public double h;


        public MGCell()
        {
            nbOfParticles = vertexCount();
        }

        public MGCell(int i)
        {
            cellId = i;
            nuclei = new Vertex[2];
            nuclei[0] = new Vertex();
            centre = new Vector(0, 0, 0);
            
            nbOfParticles = vertexCount();
            innerRadius = new Vector(1, 1, 1);
            neighbours = new List<int>();

            sigma = new int[nbOfParticles];
            spins = new float[nbOfParticles];
            targetVertices = new Vector[nbOfParticles];
            subset = new List<int>();

            externalEdges = new EdgeSet();
            nucleusEdges = new EdgeSet();

            nuclei[0].id = nbOfParticles;
            nuclei[0].cellId = cellId;

            axis = new Vector[3];
            for (int j = 0; j < nbOfParticles; j++)
            {
                vertices[j].id = j;
                vertices[j].cellId = i;
                Edge edge = new Edge(vertices[j], nuclei[0]);
                nucleusEdges.add(edge);
                targetVertices[j] = vertices[j].Clone().v;

                int k = new Random().Next(2);
                spins[j] = (k == 0) ? MGModel.delta : -MGModel.delta;
                sigma[j] = j;
                vertices[j].externalNeighbours = new List<int[]>();
                vertices[j].globalForces = new Vector();
                subset.Add(j);
            }

            
            nuclei[0].v = ComputeCentreFromMesh();
            SetEdgeELengths();
            ResetCell();
        }

        public MGCell(int i, Mesh m)
        {
            cellId = i;
            nuclei = new Vertex[2];
            nuclei[0] = new Vertex();
            centre = new Vector(0, 0, 0);
            
            ShapeCell(m);
            nbOfParticles = vertexCount();
            innerRadius = m.innerRadius;
            neighbours = new List<int>();

            sigma = new int[nbOfParticles];
            spins = new float[nbOfParticles];
            targetVertices = new Vector[nbOfParticles];

            externalEdges = new EdgeSet();
            nucleusEdges = new EdgeSet();
            subset = new List<int>();

            nuclei[0].id = nbOfParticles;
            nuclei[0].cellId = cellId;
            axis = new Vector[3];
            for (int j = 0; j < nbOfParticles; j++)
            {
                vertices[j].id = j;
                vertices[j].cellId = i;
                Edge edge = new Edge(vertices[j], nuclei[0]);
                nucleusEdges.add(edge);
                targetVertices[j] = vertices[j].Clone().v;

                int k = new Random().Next(2);
                spins[j] = (k == 0) ? MGModel.delta : -MGModel.delta;
                sigma[j] = j;
                vertices[j].externalNeighbours = new List<int[]>();
                vertices[j].globalForces = new Vector();
                subset.Add(j);
            }

            nuclei[0].v = ComputeCentreFromMesh();
            SetEdgeELengths();
            ResetCell();
        }

        public MGCell(int i, int tissueId, Mesh m)
        {
            cellId = i;
            nuclei = new Vertex[2];
            nuclei[0] = new Vertex();
            centre = new Vector(0, 0, 0);

            ShapeCell(m);
            nbOfParticles = vertexCount();
            innerRadius = m.innerRadius;
            neighbours = new List<int>();

            sigma = new int[nbOfParticles];
            spins = new float[nbOfParticles];
            targetVertices = new Vector[nbOfParticles];

            externalEdges = new EdgeSet();
            nucleusEdges = new EdgeSet();
            subset = new List<int>();

            nuclei[0].id = nbOfParticles;
            nuclei[0].cellId = cellId;
            axis = new Vector[3];
            for (int j = 0; j < nbOfParticles; j++)
            {
                vertices[j].id = j;
                vertices[j].cellId = i;
                vertices[j].tissueId = tissueId;

                Edge edge = new Edge(vertices[j], nuclei[0]);
                nucleusEdges.add(edge);
                targetVertices[j] = vertices[j].Clone().v;

                int k = new Random().Next(2);
                spins[j] = (k == 0) ? MGModel.delta : -MGModel.delta;
                sigma[j] = j;
                vertices[j].externalNeighbours = new List<int[]>();
                vertices[j].globalForces = new Vector();
                subset.Add(j);
            }

            nuclei[0].v = ComputeCentreFromMesh();
            SetEdgeELengths();
            ResetCell();
        }

        public void PositionCell(Vector c) {
            centre.Copy(c);
            nuclei[0].v.Copy(c);

            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].translate(c);
            }
        }

        public void Translate(Vector u)
        {
            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].translate(u);
            }
            nuclei[0].translate(u);
            centre += u;
        }
        
        public void ShapeCell(Mesh m)
        {
            vertices = new VertexSet();
            edges = new EdgeSet();
            faces = new FaceSet();

            copy(m);
        }
        
        public void ComputeForces()
        {
            ResetCell();
            
            for (int i = 0; i < edgeCount(); i++)
            {
                //edges[i].force = MGModel.Force(edges[i].length(), MGModel.u0, edges[i].l0);
                edges[i].force = MGModel.Force(edges[i].length(), MGModel.J[tissueId, tissueId], edges[i].l0);
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                //nucleusEdges[i].force = MGModel.Force(nucleusEdges[i].length(), MGModel.u0, nucleusEdges[i].l0);
                nucleusEdges[i].force = MGModel.Force(nucleusEdges[i].length(), MGModel.J[tissueId, tissueId], nucleusEdges[i].l0);
                //Console.WriteLine(cellId + ", " + i + ", " + nucleusEdges[i].force + "; " + nucleusEdges[i].l0 + "; " + nucleusEdges[i].length());
            }
            
            for (int i = 0; i < edgeCount(); i++)
            {
                edges[i].ends[0].internalForces += edges[i].force * edges[i].UnitVector();
                edges[i].ends[1].internalForces -= edges[i].force * edges[i].UnitVector();
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                nucleusEdges[i].ends[0].nucleusForce0 += nucleusEdges[i].force * nucleusEdges[i].UnitVector();
                nucleusEdges[i].ends[1].force -= nucleusEdges[i].force * nucleusEdges[i].UnitVector();
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                vertices[i].force = vertices[i].internalForces + vertices[i].nucleusForce0;
            }
        }

        public void ComputeExternalForces()
        {
            if (MGModel.elasticExternalSpring)
            {
                for (int i = 0; i < externalEdges.getCount(); i++)
                {
                    externalEdges[i].force = MGModel.Force(externalEdges[i].length(), MGModel.J[externalEdges[i].ends[0].tissueId, externalEdges[i].ends[1].tissueId], MGModel.DCol);
                }

                for (int i = 0; i < externalEdges.getCount(); i++)
                {
                    externalEdges[i].ends[0].externalForces += externalEdges[i].force * externalEdges[i].UnitVector();
                    externalEdges[i].ends[1].externalForces -= externalEdges[i].force * externalEdges[i].UnitVector();
                }
            }
            
            for (int i = 0; i < externalEdges.getCount(); i++)
            {
                //Console.WriteLine(cellId + ", " + i + ", " + externalEdges[i].ends[1].nucleusForce0);
                if (!externalEdges[i].ends[1].nullForces)
                    //externalEdges[i].ends[0].externalForces += externalEdges[i].ends[1].internalForces + externalEdges[i].ends[1].nucleusForce0 + externalEdges[i].ends[1].globalForces;
                    externalEdges[i].ends[1].externalForces += externalEdges[i].gamma * (externalEdges[i].ends[0].internalForces + externalEdges[i].ends[0].nucleusForce0 + externalEdges[i].ends[0].globalForces);
                //else
                //    externalEdges[i].ends[0].externalForces += externalEdges[i].ends[1].internalForces + externalEdges[i].ends[1].nucleusForce0;
                //if(cellId>=50)
                    //Console.WriteLine(externalEdges[i].ends[1].externalForces);
            }

            
            for (int i = 0; i < nbOfParticles; i++)
            {
                if (!vertices[i].nullForces) {
                    vertices[i].force = vertices[i].internalForces + vertices[i].nucleusForce0;

                    /*
                    if (vertices[i].force.norm() > MGModel.maxForce)
                        vertices[i].force *= MGModel.maxForce / vertices[i].force.norm();

                    if (vertices[i].externalForces.norm() > MGModel.maxForce)
                        vertices[i].externalForces *= MGModel.maxForce / vertices[i].externalForces.norm();
                    //*/

                    //Console.WriteLine(cellId + ", " + i + ", " + vertices[i].externalForces);
                    vertices[i].force += vertices[i].externalForces + vertices[i].globalForces;
                    
                }
                else {
                    vertices[i].force = new Vector();
                }
                //Console.WriteLine(cellId + ", " + i + ", " + vertices[i].internalForces + "; " + vertices[i].nucleusForce0 + "; " + vertices[i].externalForces);
                //Console.WriteLine(cellId + ", " + i + ", " + vertices[i].internalForces + "; " + vertices[i].nucleusForce0);
                //Console.WriteLine(cellId + ", " + i + ", " + vertices[i].force + "; " + vertices[i].force.norm());
            }
        }
        public double metric;
        public void Dynamise()
        {
            for(int i=0; i< nbOfParticles; i++) {
                vertices[i].force = vertices[i].internalForces + vertices[i].nucleusForce0;
                vertices[i].force += vertices[i].externalForces + vertices[i].globalForces;

                vertices[i].Move();
            }
            nuclei[0].Move();
            cycleTime++;
        }

        public void CellCycle(int appliedForces, Vector normal, bool staticShape)
        {
            //Compute forces applied on cell
            //ChangeShape();
            //Console.WriteLine("Cycle");
            ComputeCentreFromMesh();
            for (int j = 0; j < vertexCount(); j++)
            {
                vertices[j].externalForces = new Vector();
                vertices[j].globalForces = new Vector();
            }
            ComputeForces();

            //Cell Growth
            if (inGrowMode)
            {
                Dynamise();
            }

            //Check if cell is ready for division
            if (cycleTime == MGModel.cellCyclePeriod && !hasDivided)
            {
                inDivisionMode = true;
                inGrowMode = false;
            }

            //Cell Division
            if (inDivisionMode && Simulator.cellPopulation.populationSize < Simulator.cellPopulation.maxPopulationSize)
            {
                //Console.WriteLine("Potentially Dividing ...");
                double r = new Random().NextDouble();
                if(r < MGModel.divisionRate)
                {
                    Console.WriteLine(" Cell" + cellId + " Dividing ...");
                    Mitosis(normal, staticShape);
                    inDivisionMode = false;
                    MGModel.searchNeighbours = true;
                }
            }

            //Initialize new cells
            if (isNewBorn && hasDivided)
            {
                inGrowMode = true;
                isNewBorn = false;
                hasDivided = false;
            }
        }

        public void KleinCellCycle(int appliedForces, Vector normal, bool staticShape) {

            //Console.WriteLine(cycleTime);
            //Check if cell is ready for division
            if (cycleTime==MGModel.cellCyclePeriod && !hasDivided)
            {
                inDivisionMode = true;
                inGrowMode = false;
                //cycleTime = 0;
            }
            //Console.WriteLine("Tissue: " + tissueName + ", cycleTime: " + cycleTime);
            //Cell Division
            //if (inDivisionMode && Simulator.cellPopulation.populationSize < Simulator.cellPopulation.maxPopulationSize)
            if (inDivisionMode && Simulator.cellPopulation.tissues[tissueId].populationSize < Simulator.cellPopulation.tissues[tissueId].maxPopulationSize )
            {
                double r = new Random().NextDouble();
                //if (r <= MGModel.divisionRate)
                if (r <= divisionRate)
                {
                    //
                    Console.WriteLine(" Cell" + cellId + " Dividing ...");
                    Mitosis(normal, staticShape);
                    inDivisionMode = false;
                    MGModel.searchNeighbours = true;
                    MGModel.nextNeighboursSearchFrame = Simulator.frame;
                }
            }

            //Initialize new cells
            if (isNewBorn && hasDivided)
            {
                inGrowMode = true;
                isNewBorn = false;
                hasDivided = false;
            }
        }

        public void Mitosis(Vector normal, bool staticShape)
        {
            Plane splitPlane = CreateCutPlane2(normal);
            MeshSplitter splitter = new MeshSplitter(this, splitPlane);
            splitter.MeshInitialize();
            splitter.MeshSplit();
            Mesh lowerMesh = splitter.CreateMeshLower();
            Mesh upperMesh = splitter.CreateMeshUpper();
            
            MGCell newCell = new MGCell(Simulator.cellPopulation.populationSize, lowerMesh);
            if (newCell != null)
            {
                newCell.cellId = Simulator.cellPopulation.populationSize;
                newCell.appliedForces = appliedForces;
                newCell.Rcell = Rcell;

                newCell.inDivisionMode = false;
                newCell.inGrowMode = true;
                newCell.isNewBorn = true;
                newCell.parent = cellId;
                newCell.cyclePeriod = 200;
                newCell.cycleTime = 0;
                newCell.centre = newCell.ComputeCentreFromMesh();
                newCell.nuclei[0].v = newCell.centre;
                newCell.spins = spins;
                newCell.polarisation = polarisation;
                subset = new List<int>();
                
                for (int i = 0; i<newCell.nbOfParticles; i++)
                {
                    newCell.vertices[i].globalForces = new Vector();
                    newCell.subset.Add(i);
                }
                Simulator.cellPopulation.AddCell(tissueId, newCell);
            }

            //Booleans
            inDivisionMode = false;
            hasDivided = true;
            isNewBorn = true;
            inGrowMode = true;
            cycleTime = 0;
            //cyclePeriod = 100;

            copy(upperMesh);
            nuclei[0].v = ComputeCentreFromMesh();
            centre = ComputeCentreFromMesh();
            nucleusEdges = new EdgeSet();
            for (int i = 0; i < nbOfParticles; i++)
            {
                Edge edge = new Edge(vertices[i], nuclei[0]);
                nucleusEdges.add(edge);

                vertices[i].globalForces = new Vector();
                targetVertices[i] = vertices[i].Clone().v;
                subset.Add(i);
            }
            
            if (staticShape)
            {
                SetEdgeELengths();
                newCell.SetEdgeELengths();
            }
            else
            {
                if(polarisation == Vector.down)
                {
                    SetElengths(MGModel.modelMeshes[1]);
                    newCell.SetElengths(MGModel.modelMeshes[1]);
                }
                else
                {
                    SetElengths(MGModel.modelMeshes[0]);
                    newCell.SetElengths(MGModel.modelMeshes[0]);
                }
            }
        }

        public void ChangeShape()
        {
            double dist = 0;
            float[] minDist = new float[nbOfParticles];
            int jAfter = 0;
            int jCurrent = 0;
            Vector meanPoint = centre;
            Vector norm = Vector.zero;
            Vector tang = Vector.zero;
            Vector tP;

            Vector spinSiteCoordinatesInCurrentConfig, spinSiteCoordinatesAfterFlip;
            double energyContributionInCurrentConfig = 100f, energyContributionAfterFlip = 100f;
            Vector[] currentPos = new Vector[nbOfParticles];

            int spinSite = new Random().Next(nbOfParticles);

            for (int i = 0; i < spinSite; i++)
            {
                norm = vertices[i].GetPosition() - meanPoint;
                currentPos[i] = (nucleusEdges[i].l0 + spins[i]) * (norm / norm.norm()) + meanPoint;
            }
            for (int i = spinSite + 1; i < nbOfParticles; i++)
            {
                norm = vertices[i].GetPosition() - meanPoint;
                currentPos[i] = (nucleusEdges[i].l0 + spins[i]) * (norm / norm.norm()) + meanPoint;
            }

            norm = vertices[spinSite].GetPosition() - meanPoint;
            spinSiteCoordinatesInCurrentConfig = (nucleusEdges[spinSite].l0 + spins[spinSite]) * (norm / norm.norm()) + meanPoint;
            spinSiteCoordinatesAfterFlip = (nucleusEdges[spinSite].l0 - spins[spinSite]) * (norm / norm.norm()) + meanPoint;

            for (int j = 0; j < nbOfParticles; j++)
            {
                dist = Vector.Distance(spinSiteCoordinatesInCurrentConfig, targetVertices[j]);
                if (dist < energyContributionInCurrentConfig)
                {
                    energyContributionInCurrentConfig = dist;
                    jCurrent = j;
                }
            }

            for (int j = 0; j < nbOfParticles; j++)
            {
                dist = Vector.Distance(spinSiteCoordinatesAfterFlip, targetVertices[j]);
                if (dist < energyContributionAfterFlip)
                {
                    energyContributionAfterFlip = dist;
                    jAfter = j;
                }
            }

            double energyGap = energyContributionAfterFlip - energyContributionInCurrentConfig;

            if (energyGap < 0)
            {
                spins[spinSite] = -spins[spinSite];
                currentPos[spinSite] = spinSiteCoordinatesAfterFlip;
                tP = targetVertices[jAfter];
            }
            else
            {
                currentPos[spinSite] = spinSiteCoordinatesInCurrentConfig;
                tP = targetVertices[jCurrent];
                if (MGModel.T > 0)
                {
                    currentPos[spinSite] = spinSiteCoordinatesInCurrentConfig;
                    double p = Math.Exp(-energyGap / MGModel.T);
                    double r = new Random().NextDouble();
                    if (r <= p)
                    {
                        spins[spinSite] = -spins[spinSite];
                        currentPos[spinSite] = spinSiteCoordinatesAfterFlip;
                        tP = targetVertices[jAfter];
                    }
                }
            }

            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 += spins[i];
            }

            /*
            norm = currentPos[spinSite] - meanPoint;
            Plane plane = new Plane(currentPos[spinSite], norm);
            Vector direction = (tP - currentPos[spinSite]).normalized;
            tang = plane.PointOrthogonalProjection(direction);

            float coef = 0, max = 0; int jMax = 0;
            for (int j = 0; j < particleNeighbourhoods[spinSite].Count; j++)
            {
                coef = Vector.Dot(tang, vertices[particleNeighbourhoods[spinSite][j]].GetPosition() - currentPos[spinSite]);
                if (Math.Abs(coef) > Math.Abs(max))
                {
                    max = coef;
                    jMax = j;
                }
            }
            //if(Mathf.Abs(max) - MGModel.delta > MGModel.epsilon)
            //{
            //cellReq[spinSite][jMax] -= Mathf.Sign(max) * 0.01f;
            //Debug.Log(cellReq[spinSite][jMax]);
            //}
            */
        }

        public void Reshape(Mesh mesh)
        {
            float[,] costs = new float[vertexCount(), mesh.vertexCount()];
            for(int i=0; i<vertexCount(); i++)
            {
                for(int j=0; j<mesh.vertexCount(); j++)
                {
                    costs[i,j] = (float)Vector.Distance(vertices[i].v, mesh.vertices[j].v);
                }
            }

            int[] mapping = HungarianAlgorithm.FindAssignments(costs);

            for(int i=0; i<mapping.Length; i++)
            {
                vertices[i].v = mesh.vertices[mapping[i]].v;
            }

            //Console.WriteLine(mapping.Length);
        }

        public void ReshapeMembrane(List<int> particlesList, float scale)
        {
            for (int j = 0; j < edgeCount(); j++)
            {
                if (particlesList.Contains(edges[j].ends[0].pos) 
                        && particlesList.Contains(edges[j].ends[1].pos))
                {
                    edges[j].l0 *= scale;
                }
            }
        }

        public void SetElengths(Mesh cell)
        {
            Vector cellCentre = cell.ComputeCentreFromMesh();
            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 = Vector.Distance(cellCentre, cell.vertices[i].v);
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = cell.edges[j].length();
            }
        }

        private Plane CreateCutPlane(Vector normal)
        {
            GetElongationAxis();
            int min = (int)elongationAxis[0];
            int max = (int)elongationAxis[1];
            //Console.WriteLine(Vector.Distance(vertices[min].v, vertices[max].v));
            Vector centre = Vector.Lerp(vertices[min].GetPosition(), vertices[max].GetPosition(), .5f);

            potentialNucleiPos = new Vector[2];
            potentialNucleiPos[0] = Vector.Lerp(vertices[max].GetPosition(), centre, .5f);
            potentialNucleiPos[1] = Vector.Lerp(vertices[min].GetPosition(), centre, .5f);

            if (normal == Vector.zero)
            {
                normal = (vertices[max].GetPosition() - vertices[min].GetPosition());
                normal.normalize();
            }

            Vector newNormal = (vertices[max].GetPosition() - vertices[min].GetPosition()) ^ (ComputeCentreFromMesh() - vertices[max].GetPosition());

            return new Plane(centre, normal);
            //return new Plane(centre, newNormal);
        }

        private Plane CreateCutPlane2(Vector normal)
        {
            //centre = ComputeCentreFromMesh();
            if (normal == Vector.zero)
            {
                FindAxis();
                //normal = axis[1];
                
                Vector up = new Vector(0, 1, 0);
                double zero = Math.Abs(axis[0] * up);
                double one = Math.Abs(axis[1] * up);
                double two = Math.Abs(axis[2] * up);

                if (zero >= one && zero >= two)
                    normal = axis[0];
                if (one >= zero && one >= two)
                    normal = axis[1];
                if (two >= zero && two >= one)
                    normal = axis[2];

                normal = axis[0];
            }

            return new Plane(centre, normal);
        }

        public void FindAxis()
        {
            //centre = ComputeCentreFromMesh();
            double nCovXY = 0;
            double nCovXZ = 0;
            double nCovYZ = 0;
            double nVarX = 0;
            double nVarY = 0;
            double nVarZ = 0;

            for (int i=0; i<nbOfParticles; i++)
            {
                nCovXY += (vertices[i].v.x - centre.x) * (vertices[i].v.y - centre.y);
                nCovXZ += (vertices[i].v.x - centre.x) * (vertices[i].v.z - centre.z);
                nCovYZ += (vertices[i].v.y - centre.y) * (vertices[i].v.z - centre.z);

                nVarX += (vertices[i].v.x - centre.x) * (vertices[i].v.x - centre.x);
                nVarY += (vertices[i].v.y - centre.y) * (vertices[i].v.y - centre.y);
                nVarZ += (vertices[i].v.z - centre.z) * (vertices[i].v.z - centre.z);
            }

            Matrix3x3 M = new Matrix3x3(nVarX, nCovXY, nCovXZ,
                                        nCovXY, nVarY, nCovYZ,
                                        nCovXZ, nCovYZ, nVarZ);

            M.sym_eigen();
            M.compute_matrix();

            double minEigen = Math.Abs(M.eigenValues.x);
            int minEigenIndex = 0;
            if(Math.Abs(M.eigenValues.y) < minEigen)
            {
                minEigen = M.eigenValues.y;
                minEigenIndex = 1;
            }
            if (Math.Abs(M.eigenValues.z) < minEigen)
            {
                minEigen = M.eigenValues.z;
                minEigenIndex = 2;
            }

            double maxEigen = Math.Abs(M.eigenValues.x);
            int maxEigenIndex = 0;
            if (Math.Abs(M.eigenValues.y) > maxEigen)
            {
                maxEigen = M.eigenValues.y;
                maxEigenIndex = 1;
            }
            if (Math.Abs(M.eigenValues.z) > maxEigen)
            {
                maxEigen = M.eigenValues.z;
                maxEigenIndex = 2;
            }

            int otherEigenIndex = 3 - minEigenIndex - maxEigenIndex;

            axis[0] = new Vector(M.R[0, maxEigenIndex],
                                M.R[1, maxEigenIndex],
                                M.R[2, maxEigenIndex]);

            axis[1] = new Vector(M.R[0, otherEigenIndex],
                                M.R[1, otherEigenIndex],
                                M.R[2, otherEigenIndex]);

            axis[2] = new Vector(M.R[0, minEigenIndex],
                                M.R[1, minEigenIndex],
                                M.R[2, minEigenIndex]);
        }
        
        public int FindTopParticle2(Vector polarisation)
        {
            //FindAxis();

            double[] proj = new double[nbOfParticles];

            double max = 0;
            //double max = (vertices[0].v - centre) * polarisation;

            //centre = ComputeCentreFromMesh();

            int maxCount = 0;
            //int maxCount = 1;

            List<int> maxIndices = new List<int>();
            //maxIndices.Add(0);

            polarisation.normalize();
            //for (int i=1; i<nbOfParticles; i++)
            for (int i = 0; i < nbOfParticles; i++)
            {
                proj[i] = (vertices[i].v - centre) * polarisation;

                //if(cellId == 1)
                    //Console.WriteLine(cellId + ", " + i + ", " + proj[i] + ", " + polarisation + ", " + (vertices[i].v - centre));

                if (proj[i] == max)
                {
                    maxCount++;
                    maxIndices.Add(i);
                }

                if (proj[i] > max)
                {
                    max = proj[i];
                    maxCount = 1;
                    maxIndices.Clear();
                    maxIndices.Add(i);
                }
            }

            double min =((vertices[maxIndices[0]].v - centre) ^ polarisation).norm();
            int minIndex = maxIndices[0];
            //Console.WriteLine(cellId + ": " + maxCount);

            for (int i=0; i<maxCount; i++)
            {
                //Console.WriteLine(vertices[maxIndices[i]].pos);
                double height = ((vertices[maxIndices[i]].v - centre) ^ polarisation).norm();

                if(min > height)
                {
                    min = height;
                    minIndex = maxIndices[i];
                }
            }

            return minIndex;
        }

        public int FindTopParticle(Vector point)
        {
            double distance = 0;

            double min = 100000;
            
            int minIndex = 0;
            
            for (int i = 0; i < nbOfParticles; i++)
            {
                distance = (vertices[i].v - point).norm();
                //Console.WriteLine(cellId + ", " + i + ", " + distance);

                if (distance < min)
                {
                    min = distance;
                    minIndex = i;
                }
            }
            return minIndex;
        }

        public void GetElongationAxis()
        {
            ComputeParticleDistances();
            elongationAxis = Helper.DetermineElongationAxis(particlesDistances, nbOfParticles);
        }

        public void ComputeParticleDistances()
        {
            particlesDistances = new float[nbOfParticles, nbOfParticles];
            for (int i = 0; i < nbOfParticles; i++)
            {
                int j = 0;
                while (j++ < i)
                {
                    particlesDistances[i, j] = (float)Vector.Distance(vertices[j].GetPosition(), vertices[i].GetPosition());
                    particlesDistances[j, i] = particlesDistances[i, j];
                }
            }
        }
        
        public void ResetCell()
        {
            //externalEdges.clear();
            nuclei[0].force = new Vector();
            Parallel.For(0, nbOfParticles, i =>
            {
                vertices[i].internalForces = new Vector();
                //vertices[i].externalForces = new Vector();
                vertices[i].nucleusForce0 = new Vector();
                vertices[i].force = new Vector();
            });
        }
        
        public Vector GetPosition() {
            return centre;
            //return nuclei[0].GetPosition();
            //return ComputeCentreFromMesh();
        }
        
        public virtual void ApicalConstriction(double d){}
        public virtual void PlanarPolarisedConstriction(double d){}
        public virtual void ChangeToRBCShape(double R, double a0, double a1, double a2){}

        public void SetEdgeELengths()
        {
            for(int i=0; i < edgeCount(); i++)
            {
                edges[i].l0 = edges[i].length();
                //Console.WriteLine(edges[i].l0);
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                nucleusEdges[i].l0 = nucleusEdges[i].length();
                //Console.WriteLine(nucleusEdges[i].l0);
            }
            for (int i = 0; i < externalEdges.getCount(); i++)
            {
                externalEdges[i].l0 = MGModel.DCol;
            }
        }

        public void ScaleEdgeELengths(float scale)
        {
            for (int i = 0; i < edgeCount(); i++)
            {
                edges[i].l0 *= scale;
                //Console.WriteLine(edges[i].l0);
            }
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                nucleusEdges[i].l0 *= scale;
                //Console.WriteLine(nucleusEdges[i].l0);
            }
        }

        public void ComputeTargetELengths(){
            for (int i = 0; i < nbOfParticles; i++)
            {
                nucleusEdges[i].l0 = targetVertices[i].norm();
            }

            for (int j = 0; j < edgeCount(); j++)
            {
                edges[j].l0 = Vector.Distance(targetVertices[edges[j].ends[0].pos], targetVertices[edges[j].ends[1].pos]);
            }
        }

        public float ElasticEnergyMembraneRays()
        {
            double E=0;
            for(int i=0; i<edgeCount(); i++)
            {
                E += (edges[i].l0 - edges[i].length())* (edges[i].l0 - edges[i].length());
            }
            return (float)E;
        }

        public float ElasticEnergyNucleusRays()
        {
            double E = 0;
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                E += (nucleusEdges[i].l0 - nucleusEdges[i].length()) * (nucleusEdges[i].l0 - nucleusEdges[i].length());
            }
            return (float)E;
        }

        public float MorseEnergyMembraneRays()
        {
            double E = 0;
            for (int i = 0; i < edgeCount(); i++)
            {
                double exponent = MGModel.rho * (edges[i].l0 - edges[i].length());
                E += MGModel.J[tissueId, tissueId] * (Math.Exp(2 * exponent) - MGModel.alpha * Math.Exp(exponent));
            }
            return (float)E;
        }

        public float MorseEnergyNucleusRays()
        {
            double E = 0;
            for (int i = 0; i < nucleusEdges.getCount(); i++)
            {
                double exponent = MGModel.rho * (nucleusEdges[i].l0 - nucleusEdges[i].length());
                E += MGModel.J[tissueId, tissueId] * (Math.Exp(2 * exponent) - MGModel.alpha * Math.Exp(exponent));
            }
            return (float)E;
        }

        public void Randomize()
        {
            Randomizer.Randomize<int>(sigma);
        }
    }
}