// Fill out your copyright notice in the Description page of Project Settings.

#include "pltg/DelaunayTriangulation.h"
#include "DrawDebugHelpers.h"
#include "Engine/Engine.h"
#include "Engine/World.h"

// Sets default values
ADelaunayTriangulation::ADelaunayTriangulation()
{
    PrimaryActorTick.bCanEverTick = true;
}

// Called when the game starts or when spawned
void ADelaunayTriangulation::BeginPlay()
{
    Super::BeginPlay();
}

// Called every frame
void ADelaunayTriangulation::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);

    // Draw debug visualization each frame
    if (bDrawTriangulation)
    {
        DrawDebugTriangulation();
    }
}

bool FDTriangle::operator==(const FDTriangle& Other) const
{
    // Create sorted arrays of vertices for comparison
    TArray<FVector2D> ThisVertices = { Vertex1, Vertex2, Vertex3 };
    TArray<FVector2D> OtherVertices = { Other.Vertex1, Other.Vertex2, Other.Vertex3 };

    // Sort vertices by X then Y
    ThisVertices.Sort([](const FVector2D& A, const FVector2D& B) {
        if (A.X == B.X) return A.Y < B.Y;
        return A.X < B.X;
        });

    OtherVertices.Sort([](const FVector2D& A, const FVector2D& B) {
        if (A.X == B.X) return A.Y < B.Y;
        return A.X < B.X;
        });

    // Compare sorted vertices
    return ThisVertices[0] == OtherVertices[0] &&
        ThisVertices[1] == OtherVertices[1] &&
        ThisVertices[2] == OtherVertices[2];
}

bool FDEdge::operator==(const FDEdge& Other) const
{
    // Edges are equal regardless of direction
    return (Start == Other.Start && End == Other.End) ||
        (Start == Other.End && End == Other.Start);
}

void ADelaunayTriangulation::GenerateTriangulation(const TArray<FVector2D>& Points)
{
    if (Points.Num() < 3)
    {
        UE_LOG(LogTemp, Warning, TEXT("Need at least 3 points for triangulation"));
        return;
    }

    // Clear previous triangulation
    Edges.Empty();
    Triangles.Empty();

    // Generate triangulation using Bowyer-Watson algorithm
    Triangles = BowyerWatson(Points);

    // Extract unique edges from triangles
    Edges = GetUniqueEdges(Triangles);

    UE_LOG(LogTemp, Log, TEXT("Generated triangulation with %d triangles and %d edges"), Triangles.Num(), Edges.Num());
}

TArray<FDTriangle> ADelaunayTriangulation::BowyerWatson(const TArray<FVector2D>& Points)
{
    TArray<FDTriangle> ResultTriangles; // Изменили имя переменной

    // Find bounding box of points
    FVector2D MinPoint(FLT_MAX, FLT_MAX);
    FVector2D MaxPoint(-FLT_MAX, -FLT_MAX);

    for (const FVector2D& Point : Points)
    {
        MinPoint.X = FMath::Min(MinPoint.X, Point.X);
        MinPoint.Y = FMath::Min(MinPoint.Y, Point.Y);
        MaxPoint.X = FMath::Max(MaxPoint.X, Point.X);
        MaxPoint.Y = FMath::Max(MaxPoint.Y, Point.Y);
    }

    // Calculate dimensions
    float DX = MaxPoint.X - MinPoint.X;
    float DY = MaxPoint.Y - MinPoint.Y;
    float DeltaMax = FMath::Max(DX, DY);
    FVector2D MidPoint = (MinPoint + MaxPoint) * 0.5f;

    // Create super triangle that contains all points
    FDTriangle SuperTriangle(
        FVector2D(MidPoint.X - 20 * DeltaMax, MidPoint.Y - 10 * DeltaMax),
        FVector2D(MidPoint.X, MidPoint.Y + 20 * DeltaMax),
        FVector2D(MidPoint.X + 20 * DeltaMax, MidPoint.Y - 10 * DeltaMax)
    );

    ResultTriangles.Add(SuperTriangle); // Используем новое имя

    // Add points one by one
    for (const FVector2D& Point : Points)
    {
        TArray<FDTriangle> BadTriangles;

        // Find all triangles that are no longer valid due to the new point
        for (const FDTriangle& Triangle : ResultTriangles) // Используем новое имя
        {
            if (IsPointInCircumcircle(Triangle, Point))
            {
                BadTriangles.Add(Triangle);
            }
        }

        // Find the boundary of the polygonal hole
        TArray<FDEdge> Polygon;

        for (const FDTriangle& Triangle : BadTriangles)
        {
            // Get edges of this triangle
            TArray<FDEdge> TriangleEdges = { // Изменили имя переменной
                FDEdge(Triangle.Vertex1, Triangle.Vertex2),
                FDEdge(Triangle.Vertex1, Triangle.Vertex3),
                FDEdge(Triangle.Vertex2, Triangle.Vertex3)
            };

            // Check each edge
            for (const FDEdge& Edge : TriangleEdges) // Используем новое имя
            {
                bool IsShared = false;

                // Check if this edge is shared with any other bad triangle
                for (const FDTriangle& OtherTriangle : BadTriangles)
                {
                    if (Triangle == OtherTriangle) continue;

                    TArray<FDEdge> OtherEdges = {
                        FDEdge(OtherTriangle.Vertex1, OtherTriangle.Vertex2),
                        FDEdge(OtherTriangle.Vertex1, OtherTriangle.Vertex3),
                        FDEdge(OtherTriangle.Vertex2, OtherTriangle.Vertex3)
                    };

                    for (const FDEdge& OtherEdge : OtherEdges)
                    {
                        if (Edge == OtherEdge)
                        {
                            IsShared = true;
                            break;
                        }
                    }

                    if (IsShared) break;
                }

                // If edge is not shared, add it to the polygon
                if (!IsShared)
                {
                    Polygon.Add(Edge);
                }
            }
        }

        // Remove bad triangles from triangulation
        for (const FDTriangle& BadTriangle : BadTriangles)
        {
            ResultTriangles.Remove(BadTriangle); // Используем новое имя
        }

        // Create new triangles from the point to each edge of the polygon
        for (const FDEdge& Edge : Polygon)
        {
            ResultTriangles.Add(FDTriangle(Edge.Start, Edge.End, Point)); // Используем новое имя
        }
    }

    // Remove triangles that contain vertices from the super triangle
    TArray<FVector2D> SuperTriangleVertices = {
        SuperTriangle.Vertex1,
        SuperTriangle.Vertex2,
        SuperTriangle.Vertex3
    };

    ResultTriangles.RemoveAll([&SuperTriangleVertices](const FDTriangle& Triangle) { // Используем новое имя
        for (const FVector2D& Vertex : SuperTriangleVertices)
        {
            if (Triangle.Vertex1 == Vertex || Triangle.Vertex2 == Vertex || Triangle.Vertex3 == Vertex)
            {
                return true;
            }
        }
        return false;
        });

    return ResultTriangles; // Используем новое имя
}

bool ADelaunayTriangulation::IsPointInCircumcircle(const FDTriangle& Triangle, const FVector2D& Point)
{
    FVector2D A = Triangle.Vertex1;
    FVector2D B = Triangle.Vertex2;
    FVector2D C = Triangle.Vertex3;

    // Calculate determinant
    float D = (A.X * (B.Y - C.Y) + B.X * (C.Y - A.Y) + C.X * (A.Y - B.Y)) * 2;
    if (FMath::Abs(D) < SMALL_NUMBER) return false;

    // Calculate center of circumcircle
    float ALengthSq = A.X * A.X + A.Y * A.Y;
    float BLengthSq = B.X * B.X + B.Y * B.Y;
    float CLengthSq = C.X * C.X + C.Y * C.Y;

    float Ux = (ALengthSq * (B.Y - C.Y) + BLengthSq * (C.Y - A.Y) + CLengthSq * (A.Y - B.Y)) / D;
    float Uy = (ALengthSq * (C.X - B.X) + BLengthSq * (A.X - C.X) + CLengthSq * (B.X - A.X)) / D;

    FVector2D Center(Ux, Uy);

    // Calculate radius squared
    float RadiusSq = FVector2D::DistSquared(Center, A);

    // Calculate distance from point to center
    float DistSq = FVector2D::DistSquared(Center, Point);

    // Check if point is inside circumcircle
    return DistSq <= RadiusSq;
}

TArray<FDEdge> ADelaunayTriangulation::GetUniqueEdges(const TArray<FDTriangle>& InTriangles) // Добавили параметр
{
    TArray<FDEdge> UniqueEdges;

    for (const FDTriangle& Triangle : InTriangles) // Используем параметр
    {
        TArray<FDEdge> TriangleEdges = {
            FDEdge(Triangle.Vertex1, Triangle.Vertex2),
            FDEdge(Triangle.Vertex1, Triangle.Vertex3),
            FDEdge(Triangle.Vertex2, Triangle.Vertex3)
        };

        for (const FDEdge& Edge : TriangleEdges)
        {
            bool bFound = false;

            // Check if edge already exists (in any direction)
            for (const FDEdge& ExistingEdge : UniqueEdges)
            {
                if (Edge == ExistingEdge)
                {
                    bFound = true;
                    break;
                }
            }

            if (!bFound)
            {
                UniqueEdges.Add(Edge);
            }
        }
    }

    return UniqueEdges;
}

void ADelaunayTriangulation::DrawDebugTriangulation()
{
    UWorld* World = GetWorld();
    if (!World) return;

    // Draw all edges
    for (const FDEdge& Edge : Edges)
    {
        FVector Start(Edge.Start.X, Edge.Start.Y, 0);
        FVector End(Edge.End.X, Edge.End.Y, 0);

        DrawDebugLine(World, Start, End, EdgeColor, false, -1, 0, EdgeThickness);
    }

    // Draw points at vertices
    for (const FDEdge& Edge : Edges)
    {
        FVector StartPos(Edge.Start.X, Edge.Start.Y, 0);
        FVector EndPos(Edge.End.X, Edge.End.Y, 0);

        DrawDebugSphere(World, StartPos, DebugSphereSize, 12, EdgeColor, false, -1, 0);
        DrawDebugSphere(World, EndPos, DebugSphereSize, 12, EdgeColor, false, -1, 0);
    }
}

double ADelaunayTriangulation::Determinant(const FVector2D& A, const FVector2D& B, const FVector2D& C)
{
    return (B.X - A.X) * (C.Y - A.Y) - (B.Y - A.Y) * (C.X - A.X);
}