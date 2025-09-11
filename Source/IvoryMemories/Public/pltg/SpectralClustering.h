#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "pltg/DelaunayTriangulation.h" // Добавляем включение для FDEdge
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

    // Добавляем сеттеры
    UFUNCTION(BlueprintCallable, Category = "Spectral Clustering")
    void SetbDrawClusters(bool bDraw) { bDrawClusters = bDraw; }

    UFUNCTION(BlueprintCallable, Category = "Spectral Clustering")
    void SetDebugSphereRadius(float Radius) { DebugSphereRadius = Radius; }

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering")
    int32 NumberOfClusters = 3;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering")
    TArray<FLinearColor> ClusterColors;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering|Debug")
    bool bDrawClusters = true;

    UPROPERTY(EditAnywhere, Category = "Spectral Clustering|Debug")
    float DebugSphereRadius = 100.0f;

private:
    TArray<int32> ClusterAssignments;

    TArray<TArray<float>> BuildAdjacencyMatrix(const TArray<FVector2D>& Vertices, const TArray<FDEdge>& Edges);
    TArray<TArray<float>> ComputeDegreeMatrix(const TArray<TArray<float>>& AdjacencyMatrix);
    TArray<TArray<float>> ComputeNormalizedLaplacian(const TArray<TArray<float>>& AdjacencyMatrix, const TArray<TArray<float>>& DegreeMatrix);
    TArray<TArray<float>> ComputeEigenvectors(const TArray<TArray<float>>& LaplacianMatrix, int32 k);
    TArray<int32> KMeansClustering(const TArray<TArray<float>>& Eigenvectors, int32 k, int32 MaxIterations = 100);

    // store the vertices used for clustering
    TArray<FVector2D> ClusteredVertices;
};