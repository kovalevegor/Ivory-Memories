#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "pltg/DelaunayTriangulation.h"
#include "MinimumSpanningTree.generated.h"

// Structure representing a weighted edge
USTRUCT(BlueprintType)
struct FWeightedEdge
{
    GENERATED_BODY()

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly)
    FVector2D Start;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly)
    FVector2D End;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly)
    float Weight;

    FWeightedEdge() : Start(FVector2D::ZeroVector), End(FVector2D::ZeroVector), Weight(0.0f) {}
    FWeightedEdge(FVector2D InStart, FVector2D InEnd, float InWeight)
        : Start(InStart), End(InEnd), Weight(InWeight) {
    }
};

UCLASS()
class IVORYMEMORIES_API AMinimumSpanningTree : public AActor
{
    GENERATED_BODY()

public:
    AMinimumSpanningTree();

protected:
    virtual void BeginPlay() override;

public:
    virtual void Tick(float DeltaTime) override;

    // Function to generate MST from Delaunay edges
    UFUNCTION(BlueprintCallable, Category = "Minimum Spanning Tree")
    void GenerateMST(const TArray<FDEdge>& DelaunayEdges);

    // Function to draw debug visualization
    void DrawDebugMST();

    // Edges of the MST
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Minimum Spanning Tree")
    TArray<FDEdge> MSTEdges;

    // Debug drawing properties
    UPROPERTY(EditAnywhere, Category = "Minimum Spanning Tree|Debug")
    bool bDrawMST = true;

    UPROPERTY(EditAnywhere, Category = "Minimum Spanning Tree|Debug")
    FColor MSTColor = FColor::Blue;

    UPROPERTY(EditAnywhere, Category = "Minimum Spanning Tree|Debug")
    float MSTThickness = 3.0f;

    // Probability to add random edges back to MST (for more connected graphs)
    UPROPERTY(EditAnywhere, Category = "Minimum Spanning Tree")
    float RandomEdgeProbability = 0.1f;

private:
    // Disjoint Set Union (DSU) for Kruskal's algorithm
    class DisjointSetUnion
    {
    public:
        void MakeSet(const FVector2D& Point);
        FVector2D Find(const FVector2D& Point);
        void Unite(const FVector2D& A, const FVector2D& B);

    private:
        TMap<FVector2D, FVector2D> Parent;
        TMap<FVector2D, int32> Rank;
    };

    // Helper function to calculate distance between two points
    float CalculateDistance(const FVector2D& A, const FVector2D& B);

    // Helper function to add random edges back to MST
    TArray<FDEdge> AddRandomEdges(const TArray<FDEdge>& InMSTEdges, const TArray<FDEdge>& AllEdges, float Probability);
};