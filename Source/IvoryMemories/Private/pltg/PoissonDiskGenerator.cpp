// Fill out your copyright notice in the Description page of Project Settings.

#include "pltg/PoissonDiskGenerator.h"
#include "pltg/DelaunayTriangulation.h"
#include "DrawDebugHelpers.h"
#include "Kismet/KismetMathLibrary.h"
#include "Engine/Engine.h"

//class ADelaunayTriangulation;

// Sets default values
APoissonDiskGenerator::APoissonDiskGenerator()
{
    PrimaryActorTick.bCanEverTick = false;

    // Generate points in editor
#if WITH_EDITOR
    if (!IsRunningGame())
    {
        GeneratePoints();
    }
#endif

    TriangulationClass = ADelaunayTriangulation::StaticClass();
    TriangulationActor = nullptr;
}

void APoissonDiskGenerator::BeginPlay()
{
    Super::BeginPlay();
    GeneratePoints();
}

#if WITH_EDITOR
void APoissonDiskGenerator::PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent)
{
    Super::PostEditChangeProperty(PropertyChangedEvent);
    GeneratePoints();
}
#endif

void APoissonDiskGenerator::GeneratePoints()
{
    GeneratedPoints.Empty();
    const int N = 2;
    TArray<FVector2D> ActiveList;

    float CellSize = FMath::Max(Radius / FMath::Sqrt(static_cast<float>(N)), 1.0f);
    int GridWidth = FMath::CeilToInt(Width / CellSize) + 1;
    int GridHeight = FMath::CeilToInt(Height / CellSize) + 1;

    TArray<TArray<FVector2D>> Grid;
    Grid.Init(TArray<FVector2D>(), GridWidth);
    for (int i = 0; i < GridWidth; ++i)
    {
        Grid[i].Init(FVector2D(-1, -1), GridHeight);
    }

    FVector2D InitialPoint(
        FMath::FRandRange(0, Width),
        FMath::FRandRange(0, Height)
    );

    GeneratedPoints.Add(InitialPoint);
    ActiveList.Add(InitialPoint);

    int InitialX = FMath::FloorToInt(InitialPoint.X / CellSize);
    int InitialY = FMath::FloorToInt(InitialPoint.Y / CellSize);
    Grid[InitialX][InitialY] = InitialPoint;

    while (ActiveList.Num() > 0)
    {
        int RandomIndex = FMath::RandRange(0, ActiveList.Num() - 1);
        FVector2D Point = ActiveList[RandomIndex];
        bool Found = false;

        for (int i = 0; i < K; ++i)
        {
            float Angle = FMath::FRandRange(0, 2 * PI);
            float NewRadius = FMath::FRandRange(Radius, 2 * Radius);
            FVector2D NewPoint(
                Point.X + NewRadius * FMath::Cos(Angle),
                Point.Y + NewRadius * FMath::Sin(Angle)
            );

            if (NewPoint.X < 0 || NewPoint.X >= Width ||
                NewPoint.Y < 0 || NewPoint.Y >= Height)
            {
                continue;
            }

            int XIndex = FMath::FloorToInt(NewPoint.X / CellSize);
            int YIndex = FMath::FloorToInt(NewPoint.Y / CellSize);
            bool Valid = true;

            for (int x = FMath::Max(0, XIndex - 2); x <= FMath::Min(GridWidth - 1, XIndex + 2); ++x)
            {
                for (int y = FMath::Max(0, YIndex - 2); y <= FMath::Min(GridHeight - 1, YIndex + 2); ++y)
                {
                    if (Grid[x][y] != FVector2D(-1, -1))
                    {
                        float Distance = FVector2D::Distance(NewPoint, Grid[x][y]);
                        if (Distance < Radius)
                        {
                            Valid = false;
                            break;
                        }
                    }
                }
                if (!Valid) break;
            }

            if (Valid)
            {
                GeneratedPoints.Add(NewPoint);
                ActiveList.Add(NewPoint);
                Grid[XIndex][YIndex] = NewPoint;
                Found = true;
            }
        }

        if (!Found)
        {
            ActiveList.RemoveAt(RandomIndex);
        }

        
    }

    // Create triangulation if needed
    if (TriangulationClass && !TriangulationActor && GetWorld())
    {
        FActorSpawnParameters SpawnParams;
        SpawnParams.Owner = this;
        TriangulationActor = GetWorld()->SpawnActor<ADelaunayTriangulation>(TriangulationClass, GetActorLocation(), GetActorRotation(), SpawnParams);
    }

    // Generate triangulation
    if (TriangulationActor)
    {
        TriangulationActor->GenerateTriangulation(GeneratedPoints);
    }

    DrawDebugSpheres();
}

void APoissonDiskGenerator::DrawDebugSpheres()
{
    UWorld* World = GetWorld();
    if (!World) return;

    // Use FlushPersistentDebugLines instead of FlushDebugLines
    FlushPersistentDebugLines(World);

    for (const FVector2D& Point : GeneratedPoints)
    {
        FVector Location(Point.X, Point.Y, 0);
        DrawDebugSphere(World, Location, SphereRadius, SphereSegments, SphereColor, bPersistent);
    }
}