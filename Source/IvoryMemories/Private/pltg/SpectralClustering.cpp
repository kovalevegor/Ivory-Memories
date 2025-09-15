#include "pltg/SpectralClustering.h"
#include "pltg/DelaunayTriangulation.h"
#include "DrawDebugHelpers.h"
#include "Kismet/KismetMathLibrary.h"
#include "Engine/World.h"
#include "Math/UnrealMathUtility.h"
#include "Math/RandomStream.h"
#include "Misc/DateTime.h"

ASpectralClustering::ASpectralClustering()
{
    PrimaryActorTick.bCanEverTick = true;
    bNeedToDrawClusters = false;

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

    if (bDrawClusters && bNeedToDrawClusters)
    {
        DrawDebugClusters();
        bNeedToDrawClusters = false;
    }
}

void ASpectralClustering::PerformClustering(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges)
{
    UE_LOG(LogTemp, Log, TEXT("PerformClustering called with %d vertices and %d edges (MST + random)"),
        Vertices.Num(), Edges.Num());

    if (Vertices.Num() < NumberOfClusters)
    {
        UE_LOG(LogTemp, Warning, TEXT("Not enough vertices for clustering: %d < %d"),
            Vertices.Num(), NumberOfClusters);
        return;
    }

    // Store the vertices for debug drawing
    ClusteredVertices = Vertices;
    ClusterAssignments.Empty();

    // 1. Build similarity matrix using ALL edges (MST + random)
    UE_LOG(LogTemp, Log, TEXT("Building similarity matrix using ALL edges (MST + random)"));
    TArray<TArray<float>> SimilarityMatrix = BuildSimilarityMatrixFromEdges(Vertices, Edges);

    // 2. Build degree matrix
    UE_LOG(LogTemp, Log, TEXT("Building degree matrix"));
    TArray<TArray<float>> DegreeMatrix = BuildDegreeMatrix(SimilarityMatrix);

    // 3. Build normalized Laplacian matrix
    UE_LOG(LogTemp, Log, TEXT("Building normalized Laplacian matrix"));
    TArray<TArray<float>> LaplacianMatrix = BuildNormalizedLaplacian(SimilarityMatrix, DegreeMatrix);

    // 4. Compute eigenvectors
    UE_LOG(LogTemp, Log, TEXT("Computing eigenvectors"));
    TArray<TArray<float>> Eigenvectors;
    ComputeTopEigenvectors(LaplacianMatrix, NumberOfClusters, Eigenvectors);

    // Check if eigenvectors were computed correctly
    if (Eigenvectors.Num() == 0 || Eigenvectors[0].Num() != Vertices.Num())
    {
        UE_LOG(LogTemp, Error, TEXT("Eigenvectors computation failed or dimensions mismatch"));
        return;
    }

    // 5. Transpose and normalize eigenvectors
    UE_LOG(LogTemp, Log, TEXT("Transposing and normalizing eigenvectors"));
    TArray<TArray<float>> TransposedEigenvectors;
    TransposedEigenvectors.SetNum(Vertices.Num());

    for (int32 i = 0; i < Vertices.Num(); i++)
    {
        TransposedEigenvectors[i].SetNum(NumberOfClusters);
        for (int32 j = 0; j < NumberOfClusters; j++)
        {
            TransposedEigenvectors[i][j] = Eigenvectors[j][i];
        }

        // Normalize each row
        float Norm = 0.0f;
        for (int32 j = 0; j < NumberOfClusters; j++)
        {
            Norm += FMath::Square(TransposedEigenvectors[i][j]);
        }

        Norm = FMath::Sqrt(Norm);
        if (Norm > 0)
        {
            for (int32 j = 0; j < NumberOfClusters; j++)
            {
                TransposedEigenvectors[i][j] /= Norm;
            }
        }
    }

    // 6. Perform K-means clustering with multiple runs
    UE_LOG(LogTemp, Log, TEXT("Starting K-means clustering with %d runs"), KMeansRuns);
    ClusterAssignments = KMeansClustering(TransposedEigenvectors, NumberOfClusters, KMeansRuns);

    UE_LOG(LogTemp, Log, TEXT("Spectral clustering completed with %d clusters. Vertices: %d, Assignments: %d"),
        NumberOfClusters, ClusteredVertices.Num(), ClusterAssignments.Num());

    // Debug: count points in each cluster
    TArray<int32> ClusterCounts;
    ClusterCounts.Init(0, NumberOfClusters);

    for (int32 Assignment : ClusterAssignments)
    {
        if (ClusterCounts.IsValidIndex(Assignment))
        {
            ClusterCounts[Assignment]++;
        }
    }

    for (int32 i = 0; i < NumberOfClusters; i++)
    {
        UE_LOG(LogTemp, Log, TEXT("Cluster %d: %d points"), i, ClusterCounts[i]);
    }

    // Final check for consistency
    if (ClusteredVertices.Num() != ClusterAssignments.Num())
    {
        UE_LOG(LogTemp, Error, TEXT("CRITICAL ERROR: Vertices count (%d) doesn't match assignments count (%d)"),
            ClusteredVertices.Num(), ClusterAssignments.Num());
    }

    // Set flag to draw clusters
    bNeedToDrawClusters = true;
}

TArray<TArray<float>> ASpectralClustering::BuildSimilarityMatrixFromEdges(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges)
{
    int32 N = Vertices.Num();
    TArray<TArray<float>> SimilarityMatrix;
    SimilarityMatrix.SetNum(N);

    for (int32 i = 0; i < N; i++)
    {
        SimilarityMatrix[i].SetNum(N);
        for (int32 j = 0; j < N; j++)
        {
            SimilarityMatrix[i][j] = 0.0f;
        }
    }

    // Calculate average edge length for adaptive sigma (using ALL edges)
    float TotalEdgeLength = 0.0f;
    int32 ValidEdgesCount = 0;

    for (const FDEdge& Edge : Edges)
    {
        int32 i = FindVertexIndex(Vertices, Edge.Start);
        int32 j = FindVertexIndex(Vertices, Edge.End);

        if (i != INDEX_NONE && j != INDEX_NONE && i != j)
        {
            float EdgeLength = FVector2D::Distance(Vertices[i], Vertices[j]);
            TotalEdgeLength += EdgeLength;
            ValidEdgesCount++;
        }
    }

    if (ValidEdgesCount == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("No valid edges found for similarity matrix"));
        return SimilarityMatrix;
    }

    float AvgEdgeLength = TotalEdgeLength / ValidEdgesCount;
    float AdaptiveSigma = AvgEdgeLength * SigmaCoefficient;

    UE_LOG(LogTemp, Log, TEXT("Adaptive sigma: %.2f (avg edge length: %.2f, coefficient: %.2f, edges: %d)"),
        AdaptiveSigma, AvgEdgeLength, SigmaCoefficient, ValidEdgesCount);

    float TwoSigmaSquared = 2.0f * FMath::Square(AdaptiveSigma);

    // Build similarity matrix based on ALL edges (MST + random)
    for (const FDEdge& Edge : Edges)
    {
        int32 i = FindVertexIndex(Vertices, Edge.Start);
        int32 j = FindVertexIndex(Vertices, Edge.End);

        if (i != INDEX_NONE && j != INDEX_NONE && i != j)
        {
            float DistanceSquared = FVector2D::DistSquared(Vertices[i], Vertices[j]);
            float Similarity = FMath::Exp(-DistanceSquared / TwoSigmaSquared);

            SimilarityMatrix[i][j] = Similarity;
            SimilarityMatrix[j][i] = Similarity;
        }
    }

    // Ensure self-similarity is 1.0
    for (int32 i = 0; i < N; i++)
    {
        SimilarityMatrix[i][i] = 1.0f;
    }

    return SimilarityMatrix;
}

TArray<TArray<float>> ASpectralClustering::BuildDegreeMatrix(const TArray<TArray<float>>& SimilarityMatrix)
{
    int32 N = SimilarityMatrix.Num();
    TArray<TArray<float>> DegreeMatrix;
    DegreeMatrix.SetNum(N);

    for (int32 i = 0; i < N; i++)
    {
        DegreeMatrix[i].SetNum(N);
        for (int32 j = 0; j < N; j++)
        {
            DegreeMatrix[i][j] = 0.0f;
        }

        float Sum = 0.0f;
        for (int32 j = 0; j < N; j++)
        {
            Sum += SimilarityMatrix[i][j];
        }

        DegreeMatrix[i][i] = Sum;
    }

    return DegreeMatrix;
}

TArray<TArray<float>> ASpectralClustering::BuildNormalizedLaplacian(const TArray<TArray<float>>& SimilarityMatrix, const TArray<TArray<float>>& DegreeMatrix)
{
    int32 N = SimilarityMatrix.Num();
    TArray<TArray<float>> LaplacianMatrix;
    LaplacianMatrix.SetNum(N);

    // Сначала вычислим D^(-1/2)
    TArray<float> D_sqrt_inv;
    D_sqrt_inv.SetNum(N);
    for (int32 i = 0; i < N; i++)
    {
        if (DegreeMatrix[i][i] > 0)
            D_sqrt_inv[i] = 1.0f / FMath::Sqrt(DegreeMatrix[i][i]);
        else
            D_sqrt_inv[i] = 0.0f;
    }

    // L = I - D^(-1/2) * S * D^(-1/2)
    for (int32 i = 0; i < N; i++)
    {
        LaplacianMatrix[i].SetNum(N);
        for (int32 j = 0; j < N; j++)
        {
            if (i == j)
                LaplacianMatrix[i][j] = 1.0f - D_sqrt_inv[i] * SimilarityMatrix[i][j] * D_sqrt_inv[j];
            else
                LaplacianMatrix[i][j] = -D_sqrt_inv[i] * SimilarityMatrix[i][j] * D_sqrt_inv[j];
        }
    }

    return LaplacianMatrix;
}

void ASpectralClustering::ComputeTopEigenvectors(const TArray<TArray<float>>& Matrix, int32 NumEigenvectors, TArray<TArray<float>>& Eigenvectors)
{
    int32 N = Matrix.Num();
    Eigenvectors.SetNum(NumEigenvectors);

    for (int32 i = 0; i < NumEigenvectors; i++)
    {
        // Initialize random eigenvector
        Eigenvectors[i].SetNum(N);
        for (int32 j = 0; j < N; j++)
        {
            Eigenvectors[i][j] = FMath::FRandRange(-1.0f, 1.0f);
        }

        // Orthogonalize against previous eigenvectors
        for (int32 j = 0; j < i; j++)
        {
            float DotProduct = 0.0f;
            for (int32 k = 0; k < N; k++)
            {
                DotProduct += Eigenvectors[i][k] * Eigenvectors[j][k];
            }

            for (int32 k = 0; k < N; k++)
            {
                Eigenvectors[i][k] -= DotProduct * Eigenvectors[j][k];
            }
        }

        // Power iteration
        for (int32 iter = 0; iter < PowerIterations; iter++)
        {
            // Multiply matrix by vector
            TArray<float> NewVector;
            NewVector.SetNum(N);
            for (int32 j = 0; j < N; j++) NewVector[j] = 0.0f;

            for (int32 row = 0; row < N; row++)
            {
                for (int32 col = 0; col < N; col++)
                {
                    NewVector[row] += Matrix[row][col] * Eigenvectors[i][col];
                }
            }

            // Orthogonalize against previous eigenvectors
            for (int32 j = 0; j < i; j++)
            {
                float DotProduct = 0.0f;
                for (int32 k = 0; k < N; k++)
                {
                    DotProduct += NewVector[k] * Eigenvectors[j][k];
                }

                for (int32 k = 0; k < N; k++)
                {
                    NewVector[k] -= DotProduct * Eigenvectors[j][k];
                }
            }

            // Normalize
            float Norm = 0.0f;
            for (int32 j = 0; j < N; j++)
            {
                Norm += FMath::Square(NewVector[j]);
            }

            Norm = FMath::Sqrt(Norm);
            if (Norm > 0)
            {
                for (int32 j = 0; j < N; j++)
                {
                    Eigenvectors[i][j] = NewVector[j] / Norm;
                }
            }
        }
    }
}

TArray<int32> ASpectralClustering::KMeansClustering(const TArray<TArray<float>>& Data, int32 k, int32 NumRuns)
{
    int32 N = Data.Num();
    if (N == 0) return TArray<int32>();
    int32 Dim = Data[0].Num();

    TArray<int32> BestAssignments;
    float BestVariance = FLT_MAX;

    // Run K-means multiple times and choose the best result
    for (int32 run = 0; run < NumRuns; run++)
    {
        TArray<int32> Assignments;
        Assignments.SetNum(N);
        for (int32 i = 0; i < N; i++) Assignments[i] = 0;

        // Initialize centroids randomly
        TArray<TArray<float>> Centroids;
        Centroids.SetNum(k);

        // Use FRandomStream for random numbers
        FRandomStream RandomStream(FDateTime::Now().GetTicks() + run);

        TArray<int32> Indices;
        Indices.SetNum(N);
        for (int32 i = 0; i < N; i++) Indices[i] = i;

        // Shuffle indices
        for (int32 i = N - 1; i > 0; i--)
        {
            int32 j = RandomStream.RandRange(0, i);
            int32 Temp = Indices[i];
            Indices[i] = Indices[j];
            Indices[j] = Temp;
        }

        for (int32 i = 0; i < k; i++)
        {
            Centroids[i] = Data[Indices[i]];
        }

        for (int32 iter = 0; iter < 100; iter++)
        {
            // Assign points to nearest centroid
            for (int32 i = 0; i < N; i++)
            {
                float MinDistance = FLT_MAX;
                int32 BestCluster = 0;

                for (int32 j = 0; j < k; j++)
                {
                    float Distance = 0.0f;
                    for (int32 d = 0; d < Dim; d++)
                    {
                        Distance += FMath::Square(Data[i][d] - Centroids[j][d]);
                    }

                    if (Distance < MinDistance)
                    {
                        MinDistance = Distance;
                        BestCluster = j;
                    }
                }

                Assignments[i] = BestCluster;
            }

            // Update centroids
            TArray<TArray<float>> NewCentroids;
            TArray<int32> ClusterSizes;
            NewCentroids.SetNum(k);
            ClusterSizes.SetNum(k);
            for (int32 i = 0; i < k; i++)
            {
                NewCentroids[i].SetNum(Dim);
                for (int32 d = 0; d < Dim; d++) NewCentroids[i][d] = 0.0f;
                ClusterSizes[i] = 0;
            }

            for (int32 i = 0; i < N; i++)
            {
                int32 Cluster = Assignments[i];
                ClusterSizes[Cluster]++;

                for (int32 d = 0; d < Dim; d++)
                {
                    NewCentroids[Cluster][d] += Data[i][d];
                }
            }

            for (int32 i = 0; i < k; i++)
            {
                if (ClusterSizes[i] > 0)
                {
                    for (int32 d = 0; d < Dim; d++)
                    {
                        NewCentroids[i][d] /= ClusterSizes[i];
                    }
                }
            }

            // Check for convergence
            bool Converged = true;
            for (int32 i = 0; i < k; i++)
            {
                for (int32 d = 0; d < Dim; d++)
                {
                    if (FMath::Abs(NewCentroids[i][d] - Centroids[i][d]) > 0.001f)
                    {
                        Converged = false;
                        break;
                    }
                }
                if (!Converged) break;
            }

            if (Converged) break;

            Centroids = NewCentroids;
        }

        // Calculate within-cluster variance
        float TotalVariance = 0.0f;
        for (int32 i = 0; i < N; i++)
        {
            int32 Cluster = Assignments[i];
            for (int32 d = 0; d < Dim; d++)
            {
                TotalVariance += FMath::Square(Data[i][d] - Centroids[Cluster][d]);
            }
        }

        // Keep the best result
        if (TotalVariance < BestVariance)
        {
            BestVariance = TotalVariance;
            BestAssignments = Assignments;
        }
    }

    return BestAssignments;
}

void ASpectralClustering::DrawDebugClusters()
{
    UWorld* World = GetWorld();
    if (!World) return;

    if (ClusterAssignments.Num() == 0)
    {
        return;
    }

    if (ClusteredVertices.Num() != ClusterAssignments.Num())
    {
        return;
    }

    // Use persistent debug spheres
    FlushPersistentDebugLines(World);

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

        DrawDebugSphere(World, Location, DebugSphereRadius, 12, Color, true, -1, 0);
    }
}

void ASpectralClustering::ForceRedraw()
{
    UE_LOG(LogTemp, Log, TEXT("ForceRedraw called"));
    bNeedToDrawClusters = true;
}

//int32 ASpectralClustering::FindVertexIndex(const TArray<FVector2D>& Vertices, const FVector2D& Target, float Tolerance = 0.01f)
//{
//    for (int32 i = 0; i < Vertices.Num(); i++)
//    {
//        if (FMath::IsNearlyEqual(Vertices[i].X, Target.X, Tolerance) &&
//            FMath::IsNearlyEqual(Vertices[i].Y, Target.Y, Tolerance))
//        {
//            return i;
//        }
//    }
//    return INDEX_NONE;
//}