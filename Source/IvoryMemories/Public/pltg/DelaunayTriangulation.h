// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "DelaunayTriangulation.generated.h"

// Structure representing a triangle with three vertices
USTRUCT(BlueprintType)
struct FDTriangle
{
    GENERATED_BODY()

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    FVector2D Vertex1;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    FVector2D Vertex2;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    FVector2D Vertex3;

    FDTriangle() : Vertex1(FVector2D::ZeroVector), Vertex2(FVector2D::ZeroVector), Vertex3(FVector2D::ZeroVector) {}
    FDTriangle(FVector2D V1, FVector2D V2, FVector2D V3) : Vertex1(V1), Vertex2(V2), Vertex3(V3) {}

    bool operator==(const FDTriangle& Other) const;
};

// Structure representing an edge between two points
USTRUCT(BlueprintType)
struct FDEdge
{
    GENERATED_BODY()

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    FVector2D Start;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    FVector2D End;

    FDEdge() : Start(FVector2D::ZeroVector), End(FVector2D::ZeroVector) {}
    FDEdge(FVector2D StartPoint, FVector2D EndPoint) : Start(StartPoint), End(EndPoint) {}

    bool operator==(const FDEdge& Other) const;
};

UCLASS()
class IVORYMEMORIES_API ADelaunayTriangulation : public AActor
{
    GENERATED_BODY()

public:
    // Sets default values for this actor's properties
    ADelaunayTriangulation();

protected:
    // Called when the game starts or when spawned
    virtual void BeginPlay() override;

public:
    // Called every frame
    virtual void Tick(float DeltaTime) override;

    // Function to generate triangulation from points
    UFUNCTION(BlueprintCallable, Category = "Delaunay Triangulation")
    void GenerateTriangulation(const TArray<FVector2D>& Points);

    // Function to draw debug visualization
    void DrawDebugTriangulation();

    // Edges of the triangulation
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    TArray<FDEdge> Edges;

    // Triangles of the triangulation
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Delaunay Triangulation")
    TArray<FDTriangle> Triangles;

    // Debug drawing properties
    UPROPERTY(EditAnywhere, Category = "Delaunay Triangulation|Debug")
    bool bDrawTriangulation = true;

    UPROPERTY(EditAnywhere, Category = "Delaunay Triangulation|Debug")
    FColor EdgeColor = FColor::Green;

    UPROPERTY(EditAnywhere, Category = "Delaunay Triangulation|Debug")
    float EdgeThickness = 2.0f;

    UPROPERTY(EditAnywhere, Category = "Delaunay Triangulation|Debug")
    float DebugSphereSize = 20.0f;

private:
    // Bowyer-Watson algorithm for Delaunay triangulation
    TArray<FDTriangle> BowyerWatson(const TArray<FVector2D>& Points);

    // Check if a point is inside the circumcircle of a triangle
    bool IsPointInCircumcircle(const FDTriangle& Triangle, const FVector2D& Point);

    // Get all unique edges from triangles
    TArray<FDEdge> GetUniqueEdges(const TArray<FDTriangle>& InTriangles);

    // Calculate the determinant for three points
    static double Determinant(const FVector2D& A, const FVector2D& B, const FVector2D& C);
};