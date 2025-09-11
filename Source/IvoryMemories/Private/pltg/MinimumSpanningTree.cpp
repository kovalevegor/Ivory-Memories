#include "pltg/MinimumSpanningTree.h"
#include "DrawDebugHelpers.h"
#include "Engine/World.h"
#include "Kismet/KismetMathLibrary.h"

AMinimumSpanningTree::AMinimumSpanningTree()
{
    PrimaryActorTick.bCanEverTick = true;
}

void AMinimumSpanningTree::BeginPlay()
{
    Super::BeginPlay();
}

void AMinimumSpanningTree::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);

    if (bDrawMST)
    {
        DrawDebugMST();
    }
}

void AMinimumSpanningTree::GenerateMST(const TArray<FDEdge>& DelaunayEdges)
{
    MSTEdges.Empty();

    if (DelaunayEdges.Num() < 1)
    {
        UE_LOG(LogTemp, Warning, TEXT("No edges provided for MST generation"));
        return;
    }

    // Convert FDEdge to FWeightedEdge and calculate weights
    TArray<FWeightedEdge> WeightedEdges;
    for (const FDEdge& Edge : DelaunayEdges)
    {
        float Weight = CalculateDistance(Edge.Start, Edge.End);
        WeightedEdges.Add(FWeightedEdge(Edge.Start, Edge.End, Weight));
    }

    // Sort edges by weight (Kruskal's algorithm)
    WeightedEdges.Sort([](const FWeightedEdge& A, const FWeightedEdge& B) {
        return A.Weight < B.Weight;
        });

    // Initialize DSU
    DisjointSetUnion DSU;

    // Add all points to DSU
    for (const FWeightedEdge& Edge : WeightedEdges)
    {
        DSU.MakeSet(Edge.Start);
        DSU.MakeSet(Edge.End);
    }

    // Build MST
    for (const FWeightedEdge& Edge : WeightedEdges)
    {
        if (DSU.Find(Edge.Start) != DSU.Find(Edge.End))
        {
            MSTEdges.Add(FDEdge(Edge.Start, Edge.End));
            DSU.Unite(Edge.Start, Edge.End);
        }
    }

    // Optionally add some random edges back for more connectivity
    if (RandomEdgeProbability > 0.0f)
    {
        MSTEdges = AddRandomEdges(MSTEdges, DelaunayEdges, RandomEdgeProbability);
    }

    UE_LOG(LogTemp, Log, TEXT("Generated MST with %d edges"), MSTEdges.Num());
}

void AMinimumSpanningTree::DrawDebugMST()
{
    UWorld* World = GetWorld();
    if (!World) return;

    for (const FDEdge& Edge : MSTEdges)
    {
        FVector Start(Edge.Start.X, Edge.Start.Y, 10.0f); // Slightly above the triangulation
        FVector End(Edge.End.X, Edge.End.Y, 10.0f);

        DrawDebugLine(World, Start, End, MSTColor, false, -1, 0, MSTThickness);
    }
}

float AMinimumSpanningTree::CalculateDistance(const FVector2D& A, const FVector2D& B)
{
    return FVector2D::Distance(A, B);
}

TArray<FDEdge> AMinimumSpanningTree::AddRandomEdges(const TArray<FDEdge>& InMSTEdges, const TArray<FDEdge>& AllEdges, float Probability)
{
    TArray<FDEdge> Result = InMSTEdges;

    if (AllEdges.Num() == 0 || Probability <= 0.0f)
    {
        return Result;
    }

    // Create a set of MST edges for fast lookup
    TSet<FDEdge> MSTEdgeSet;
    for (const FDEdge& Edge : InMSTEdges)
    {
        MSTEdgeSet.Add(Edge);
        // Also add the reverse edge since edges are undirected
        MSTEdgeSet.Add(FDEdge(Edge.End, Edge.Start));
    }

    // Calculate average edge length
    float TotalLength = 0.0f;
    for (const FDEdge& Edge : AllEdges)
    {
        TotalLength += CalculateDistance(Edge.Start, Edge.End);
    }
    float AverageLength = TotalLength / AllEdges.Num();
    float MaxAllowedLength = AverageLength * 1.2f; // Allow edges up to 20% longer than average

    // Random number generator
    FRandomStream RandomStream(FDateTime::Now().GetTicks());

    // Add random edges with probability
    for (const FDEdge& Edge : AllEdges)
    {
        // Skip if edge is already in MST
        if (MSTEdgeSet.Contains(Edge) || MSTEdgeSet.Contains(FDEdge(Edge.End, Edge.Start)))
        {
            continue;
        }

        // Check if we should add this edge based on probability
        if (RandomStream.FRand() < Probability)
        {
            float EdgeLength = CalculateDistance(Edge.Start, Edge.End);

            // Only add edges that aren't too long
            if (EdgeLength <= MaxAllowedLength)
            {
                Result.Add(Edge);
            }
        }
    }

    return Result;
}

// DSU Implementation
void AMinimumSpanningTree::DisjointSetUnion::MakeSet(const FVector2D& Point)
{
    Parent.Add(Point, Point);
    Rank.Add(Point, 0);
}

FVector2D AMinimumSpanningTree::DisjointSetUnion::Find(const FVector2D& Point)
{
    if (Parent[Point] != Point)
    {
        Parent[Point] = Find(Parent[Point]);
    }
    return Parent[Point];
}

void AMinimumSpanningTree::DisjointSetUnion::Unite(const FVector2D& A, const FVector2D& B)
{
    FVector2D RootA = Find(A);
    FVector2D RootB = Find(B);

    if (RootA == RootB) return;

    if (Rank[RootA] < Rank[RootB])
    {
        Parent[RootA] = RootB;
    }
    else if (Rank[RootA] > Rank[RootB])
    {
        Parent[RootB] = RootA;
    }
    else
    {
        Parent[RootB] = RootA;
        Rank[RootA]++;
    }
}