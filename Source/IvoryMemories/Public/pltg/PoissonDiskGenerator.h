// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "PoissonDiskGenerator.generated.h"

class ADelaunayTriangulation;
class AMinimumSpanningTree;
class ASpectralClustering;

UCLASS()
class IVORYMEMORIES_API APoissonDiskGenerator : public AActor
{
    GENERATED_BODY()

public:
    APoissonDiskGenerator();

protected:
    virtual void BeginPlay() override;

#if WITH_EDITOR
    virtual void PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent) override;
#endif

public:
    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    float SphereRadius = 12.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    FColor SphereColor = FColor::Red;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    int SphereSegments = 6;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    bool bPersistent = true;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    float Width = 1000.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    float Height = 1000.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    float Radius = 80.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    int K = 30;

    UPROPERTY(VisibleAnywhere, Category = "Poisson Disk Sampling")
    TArray<FVector2D> GeneratedPoints;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    TSubclassOf<ADelaunayTriangulation> TriangulationClass;

    UPROPERTY(VisibleAnywhere, Category = "Poisson Disk Sampling")
    ADelaunayTriangulation* TriangulationActor;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    TSubclassOf<AMinimumSpanningTree> MSTClass;

    UPROPERTY(VisibleAnywhere, Category = "Poisson Disk Sampling")
    AMinimumSpanningTree* MSTActor;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    TSubclassOf<ASpectralClustering> SpectralClusteringClass;

    UPROPERTY(VisibleAnywhere, Category = "Poisson Disk Sampling")
    ASpectralClustering* SpectralClusteringActor;

private:
    void GeneratePoints();
    void DrawDebugSpheres();
    void CreateSpectralClusteringActor();
};