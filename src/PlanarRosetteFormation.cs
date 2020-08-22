using System;
using System.Collections.Generic;
using MGSharp.Core.GeometricPrimitives;
using MGSharp.Core.MGCellPopulation;
using MGSharp.Core.MGModels;
using MGSharp.Core.Helpers;
using System.Diagnostics;

namespace MGSharp
{
    class PlanarRosetteFormation: Simulator
    {
        public static void Main(string[] args)
        {
            Simulator simulator = new PlanarRosetteFormation();
            Simulator.name = "PlanarRosetteFormation";
            Simulator.states = new List<PersistantVertex>();
            Simulator.numbers = new List<PersistantNumbers>();
            Simulator.finalStates = new List<PersistantVertex>();
            Simulator.finalNumbers = new List<PersistantNumbers>();
            Simulator.metrics = new List<PersistantMetrics>();

            MGModel.Rcell = 1;
            simulator.SetupSimulation();
            simulator.SetModel();
            simulator.SetInitialConditions();

            simulator.LogParameters();

            Stopwatch stopwatch = new Stopwatch();
            Console.WriteLine("Simulation in Progress ...");
            stopwatch.Start();
            simulator.Update();
            stopwatch.Stop();
            Console.WriteLine("Simulation ended succesfully in " + stopwatch.ElapsedMilliseconds + " milliseconds.");
            simulator.Log();

            Console.WriteLine("Simulation Logs have been written to the Log file.");
            Console.ReadKey();
        }
        
        public override void SetupSimulation()
        {
            nbOfSimulationSteps = 2000;
            logFrequency = 10;
            //logVTK = true;

            popSize = 8;
            popMaxSize = 64;

            Cylinder34 mesh = new Cylinder34(new Vector(.5f,1f,.5f), false);

            Tissue t1 = new Tissue(8, popMaxSize, mesh);
            List<Tissue> tissueList = new List<Tissue>() { t1 };

            nbCellTypes = tissueList.Count;
            cellPopulation = new CellPopulation(popSize, popMaxSize, tissueList, 2);

            Helper.SetNCubedInN(popMaxSize);
            cellPopulation.PositionCells(false);
        }

        public override void SetInitialConditions()
        {
            //**************** WRITE YOUR CODE BETWEEN THIS LINE ****************//
            
            //************************** AND THIS LINE **************************//
        }

        public override void Update()
        {
            while (frame < nbOfSimulationSteps)
            {
                if (frame % 10 == 0)
                {
                    Randomize();
                }
                if (frame % logFrequency == 0)
                {
                    cellPopulation.LogState();
                    cellPopulation.LogNumbers();
                    cellPopulation.LogMetrics();
                }
                cellPopulation.DynamiseSystem();
                Console.WriteLine(frame + "/" + nbOfSimulationSteps);
                frame++;
            }
            cellPopulation.LogState();
            cellPopulation.LogNumbers();
            cellPopulation.LogMetrics();
            cellPopulation.LogFinalState();
            cellPopulation.LogFinalNumbers();
        }

        public override void SetModel()
        {
            MGModel.dT = 0.02f;
            MGModel.Rcell = 1;
            MGModel.maximumNeighbourDistance = 3f * MGModel.Rcell;
            MGModel.DInt = 0.1f;
            MGModel.DCol = 0.05f;
            MGModel.rho = 1.0f;
            MGModel.alpha = 2f;
            MGModel.damping = 0.5f;
            MGModel.delta = 0.001f;
            MGModel.maxExternalNeighbours = 3;
            MGModel.elasticExternalSpring = false;
            MGModel.hexagonalNeighbourhoodForCells = false;
            MGModel.mooreNeighbourhoodForCells = false;
            MGModel.staticNeighbourhood = true;

            MGModel.J = new float[nbCellTypes + 1, nbCellTypes];
            MGModel.DInteraction = new float[nbCellTypes + 1, nbCellTypes];
            for (int i = 0; i < nbCellTypes + 1; i++)
            {
                for (int j = 0; j < nbCellTypes; j++)
                {
                    MGModel.J[i, j] = 1f;
                    MGModel.DInteraction[i, j] = MGModel.DInt;
                }
            }
        }
    }
}