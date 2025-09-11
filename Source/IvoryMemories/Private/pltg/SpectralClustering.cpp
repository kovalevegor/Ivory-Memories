#include "pltg/SpectralClustering.h"
#include "pltg/DelaunayTriangulation.h"
#include "DrawDebugHelpers.h"
#include "Kismet/KismetMathLibrary.h"
#include "Engine/World.h"

ASpectralClustering::ASpectralClustering()
{
    PrimaryActorTick.bCanEverTick = true;

    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("f75049")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("5ef6ff")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("f0b537")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("9d2bf5")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("1ded83")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("2570d4")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("fb932e")));
    ClusterColors.Add(FLinearColor::FromSRGBColor(FColor::FromHex("ffffff")));
}

void ASpectralClustering::BeginPlay()
{
    Super::BeginPlay();
}

void ASpectralClustering::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);

    if (bDrawClusters)
    {
        DrawDebugClusters();
    }
}

void ASpectralClustering::PerformClustering(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges)
{
    if (Vertices.Num() < NumberOfClusters)
    {
        UE_LOG(LogTemp, Warning, TEXT("Not enough vertices for clustering"));
        return;
    }

    // Store the vertices for debug drawing
    ClusteredVertices = Vertices;

    // Build adjacency matrix
    TArray<TArray<float>> AdjacencyMatrix = BuildAdjacencyMatrix(Vertices, Edges);
    TArray<TArray<float>> DegreeMatrix = ComputeDegreeMatrix(AdjacencyMatrix);
    TArray<TArray<float>> LaplacianMatrix = ComputeNormalizedLaplacian(AdjacencyMatrix, DegreeMatrix);
    TArray<TArray<float>> Eigenvectors = ComputeEigenvectors(LaplacianMatrix, NumberOfClusters);
    ClusterAssignments = KMeansClustering(Eigenvectors, NumberOfClusters);

    UE_LOG(LogTemp, Log, TEXT("Spectral clustering completed with %d clusters"), NumberOfClusters);
}

TArray<TArray<float>> ASpectralClustering::BuildAdjacencyMatrix(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges)
{
    int32 n = Vertices.Num();
    TArray<TArray<float>> AdjacencyMatrix;
    AdjacencyMatrix.Init(TArray<float>(), n);

    for (int32 i = 0; i < n; i++)
    {
        AdjacencyMatrix[i].Init(0.0f, n);
    }

    for (const FDEdge& Edge : Edges)
    {
        int32 i = Vertices.IndexOfByPredicate([&Edge](const FVector2D& Vertex) {
            return Vertex.Equals(Edge.Start, 0.1f);
            });

        int32 j = Vertices.IndexOfByPredicate([&Edge](const FVector2D& Vertex) {
            return Vertex.Equals(Edge.End, 0.1f);
            });

        if (i != INDEX_NONE && j != INDEX_NONE)
        {
            float Distance = FVector2D::Distance(Edge.Start, Edge.End);
            float Weight = 1.0f / FMath::Max(Distance, 0.001f);
            AdjacencyMatrix[i][j] = Weight;
            AdjacencyMatrix[j][i] = Weight;
        }
    }

    return AdjacencyMatrix;
}

TArray<TArray<float>> ASpectralClustering::ComputeDegreeMatrix(const TArray<TArray<float>>& AdjacencyMatrix)
{
    int32 n = AdjacencyMatrix.Num();
    TArray<TArray<float>> DegreeMatrix;
    DegreeMatrix.Init(TArray<float>(), n);

    for (int32 i = 0; i < n; i++)
    {
        DegreeMatrix[i].Init(0.0f, n);
        float Sum = 0.0f;

        for (int32 j = 0; j < n; j++)
        {
            Sum += AdjacencyMatrix[i][j];
        }

        DegreeMatrix[i][i] = Sum;
    }

    return DegreeMatrix;
}

TArray<TArray<float>> ASpectralClustering::ComputeNormalizedLaplacian(const TArray<TArray<float>>& AdjacencyMatrix, const TArray<TArray<float>>& DegreeMatrix)
{
    int32 n = AdjacencyMatrix.Num();
    TArray<TArray<float>> LaplacianMatrix;
    LaplacianMatrix.Init(TArray<float>(), n);

    for (int32 i = 0; i < n; i++)
    {
        LaplacianMatrix[i].Init(0.0f, n);

        for (int32 j = 0; j < n; j++)
        {
            if (i == j)
            {
                LaplacianMatrix[i][j] = 1.0f;
            }
            else if (DegreeMatrix[i][i] > 0 && DegreeMatrix[j][j] > 0)
            {
                LaplacianMatrix[i][j] = -AdjacencyMatrix[i][j] / FMath::Sqrt(DegreeMatrix[i][i] * DegreeMatrix[j][j]);
            }
        }
    }

    return LaplacianMatrix;
}

TArray<TArray<float>> ASpectralClustering::ComputeEigenvectors(const TArray<TArray<float>>& LaplacianMatrix, int32 k)
{
    int32 n = LaplacianMatrix.Num();
    TArray<TArray<float>> Eigenvectors;
    Eigenvectors.Init(TArray<float>(), n);

    for (int32 i = 0; i < n; i++)
    {
        Eigenvectors[i].Init(0.0f, k);

        for (int32 j = 0; j < k; j++)
        {
            Eigenvectors[i][j] = FMath::FRandRange(-1.0f, 1.0f);
        }
    }

    return Eigenvectors;
}

TArray<int32> ASpectralClustering::KMeansClustering(const TArray<TArray<float>>& Eigenvectors, int32 k, int32 MaxIterations)
{
    int32 n = Eigenvectors.Num();
    TArray<int32> LocalClusterAssignments;
    LocalClusterAssignments.Init(0, n);

    TArray<TArray<float>> Centroids;
    Centroids.Init(TArray<float>(), k);

    for (int32 i = 0; i < k; i++)
    {
        int32 RandomIndex = FMath::RandRange(0, n - 1);
        Centroids[i] = Eigenvectors[RandomIndex];
    }

    for (int32 Iteration = 0; Iteration < MaxIterations; Iteration++)
    {
        for (int32 i = 0; i < n; i++)
        {
            float MinDistance = FLT_MAX;
            int32 BestCluster = 0;

            for (int32 j = 0; j < k; j++)
            {
                float Distance = 0.0f;
                for (int32 d = 0; d < Eigenvectors[i].Num(); d++)
                {
                    Distance += FMath::Pow(Eigenvectors[i][d] - Centroids[j][d], 2);
                }

                if (Distance < MinDistance)
                {
                    MinDistance = Distance;
                    BestCluster = j;
                }
            }

            LocalClusterAssignments[i] = BestCluster;
        }

        TArray<TArray<float>> NewCentroids;
        TArray<int32> ClusterSizes;
        NewCentroids.Init(TArray<float>(), k);
        ClusterSizes.Init(0, k);

        for (int32 i = 0; i < k; i++)
        {
            NewCentroids[i].Init(0.0f, Eigenvectors[0].Num());
        }

        for (int32 i = 0; i < n; i++)
        {
            int32 Cluster = LocalClusterAssignments[i];
            ClusterSizes[Cluster]++;

            for (int32 d = 0; d < Eigenvectors[i].Num(); d++)
            {
                NewCentroids[Cluster][d] += Eigenvectors[i][d];
            }
        }

        for (int32 i = 0; i < k; i++)
        {
            if (ClusterSizes[i] > 0)
            {
                for (int32 d = 0; d < NewCentroids[i].Num(); d++)
                {
                    NewCentroids[i][d] /= ClusterSizes[i];
                }
            }
        }

        Centroids = NewCentroids;
    }

    return LocalClusterAssignments;
}

void ASpectralClustering::DrawDebugClusters()
{
    UWorld* World = GetWorld();
    if (!World || ClusterAssignments.Num() == 0) return;

    // Get all vertices from the graph
    if (ClusteredVertices.Num() != ClusterAssignments.Num())
    {
        UE_LOG(LogTemp, Warning, TEXT("Mismatch between vertices and cluster assignments"));
        return;
    }

    // Draw debug spheres for each vertex with cluster color
    for (int32 i = 0; i < ClusteredVertices.Num(); i++)
    {
        FVector Location(ClusteredVertices[i].X, ClusteredVertices[i].Y, 0);
        int32 ClusterIndex = ClusterAssignments[i];

        // Ensure we have a valid color
        FColor Color = FColor::White;
        if (ClusterColors.IsValidIndex(ClusterIndex))
        {
            Color = ClusterColors[ClusterIndex].ToFColor(true);
        }

        DrawDebugSphere(World, Location, DebugSphereRadius, 12, Color, false, -1, 0);
    }
}