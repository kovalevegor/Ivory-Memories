#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "pltg/DelaunayTriangulation.h"
#include "SpectralClustering.generated.h"

UCLASS()
class IVORYMEMORIES_API ASpectralClustering : public AActor
{
    GENERATED_BODY()

public:
    ASpectralClustering();

protected:
    virtual void BeginPlay() override;

public:
    virtual void Tick(float DeltaTime) override;

    UFUNCTION(BlueprintCallable, Category = "Spectral Clustering")
    void PerformClustering(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges);

    void DrawDebugClusters();

    // Сеттеры
    UFUNCTION(BlueprintCallable, Category = "Spectral Clustering")
    void SetbDrawClusters(bool bDraw) { bDrawClusters = bDraw; }

    UFUNCTION(BlueprintCallable, Category = "Spectral Clustering")
    void SetDebugSphereRadius(float Radius) { DebugSphereRadius = Radius; }

    UFUNCTION(BlueprintCallable, Category = "Spectral Clustering")
    void ForceRedraw();

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering")
    int32 NumberOfClusters = 3;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering", meta = (ClampMin = "0.1", ClampMax = "5.0"))
    float SigmaCoefficient = 1.0f; // Коэффициент для адаптивного выбора sigma

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering")
    TArray<FLinearColor> ClusterColors;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering|Debug")
    bool bDrawClusters = true;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering|Debug")
    float DebugSphereRadius = 100.0f;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering")
    int32 KMeansRuns = 5;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering")
    int32 PowerIterations = 100;

private:
    TArray<int32> ClusterAssignments;
    TArray<FVector2D> ClusteredVertices;
    bool bNeedToDrawClusters;

    // Матричные операции
    TArray<TArray<float>> BuildSimilarityMatrixFromEdges(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges);
    TArray<TArray<float>> BuildDegreeMatrix(const TArray<TArray<float>>& SimilarityMatrix);
    TArray<TArray<float>> BuildNormalizedLaplacian(const TArray<TArray<float>>& SimilarityMatrix, const TArray<TArray<float>>& DegreeMatrix);

    // Вычисление собственных векторов методом степенной итерации
    void ComputeTopEigenvectors(const TArray<TArray<float>>& Matrix, int32 NumEigenvectors, TArray<TArray<float>>& Eigenvectors);

    // K-средних для кластеризации
    TArray<int32> KMeansClustering(const TArray<TArray<float>>& Data, int32 k, int32 NumRuns);

    int32 FindVertexIndex(const TArray<FVector2D>& Vertices, const FVector2D& Target, float Tolerance = 0.01f)
    {
        for (int32 i = 0; i < Vertices.Num(); i++)
        {
            if (FMath::IsNearlyEqual(Vertices[i].X, Target.X, Tolerance) &&
                FMath::IsNearlyEqual(Vertices[i].Y, Target.Y, Tolerance))
            {
                return i;
            }
        }
        return INDEX_NONE;
    }
};