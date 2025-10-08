// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/StaticMeshComponent.h" // UStaticMeshComponent
#include "Engine/StaticMesh.h" // UStaticMesh
#include "Materials/MaterialInterface.h" // UMaterialInterface
#include "Figure.generated.h"

class ACell;

UCLASS(Abstract) // cannot be spawned directly on level
class IVORYMEMORIES_API AFigure : public AActor
{
    GENERATED_BODY()

public:
    // Sets default values for this actor's properties
    AFigure();

protected:
    // Called when the game starts or when spawned
    virtual void BeginPlay() override;

    // Can be edited in Editor
    virtual void PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent) override;

public:
    // ---PROPERTIES---
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Figure")
    FString FigureName; // unique name for each figure

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Figure")
    ACell* CurrentCell; // cell reference for bijection

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Figure")
    int32 Rank;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Figure")
    int32 MaxMoveDistance;

    // mesh properties (editable in UE Editor)
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mesh")
    UStaticMesh* FigureMeshAsset;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mesh")
    UMaterialInterface* FigureMaterial;

    // Mesh component (not editable but visible)
    UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Components")
    UStaticMeshComponent* FigureMesh;

    // ---FUNCTIONS---
    UFUNCTION(BlueprintCallable, Category = "Figure")
    virtual bool CanMoveTo(ACell* TargetCell) const; // Check if figure can move to a cell

    UFUNCTION(BlueprintCallable, Category = "Figure")
    virtual void MoveToCell(ACell* NewCell);

    UFUNCTION(BlueprintCallable, Category = "Figure")
    virtual bool CanAttack(AFigure* TargetFigure) const; // Check if figure can attack opponent figure

    UFUNCTION(BlueprintCallable, Category = "Figure")
    virtual void Attack(AFigure* TargetFigure);

    UFUNCTION(BlueprintCallable, Category = "Figure")
    void GenerateUniqueName(const FString& BaseName, int32 Index); // Generate unique name of a class child

    UFUNCTION(BlueprintCallable, Category = "Figure")
    void SetCell(ACell* NewCell);

    UFUNCTION(BlueprintCallable, CallInEditor, Category = "Mesh")
    void ApplyMeshSettings();
};