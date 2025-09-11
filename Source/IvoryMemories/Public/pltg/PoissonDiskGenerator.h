// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "PoissonDiskGenerator.generated.h"

class ADelaunayTriangulation;

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
    float SphereRadius = 100.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    FColor SphereColor = FColor::Red;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    int SphereSegments = 12;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling|Debug")
    bool bPersistent = true;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    float Width = 1000.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    float Height = 1000.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    float Radius = 50.0f;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    int K = 30;

    UPROPERTY(VisibleAnywhere, Category = "Poisson Disk Sampling")
    TArray<FVector2D> GeneratedPoints;

    UPROPERTY(EditAnywhere, Category = "Poisson Disk Sampling")
    TSubclassOf<ADelaunayTriangulation> TriangulationClass;

    UPROPERTY(VisibleAnywhere, Category = "Poisson Disk Sampling")
    ADelaunayTriangulation* TriangulationActor;

private:
    void GeneratePoints();
    void DrawDebugSpheres();
};