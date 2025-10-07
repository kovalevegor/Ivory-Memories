// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Cell.generated.h"

class AFigure; // forward declaration

UENUM(BlueprintType) // Blueprint editable
enum class ECellState : uint8
{
	Empty UMETA(DisplayName = "Empty Cell"),
	Occupied UMETA(DisplayName = "Occupied Cell"),
	Inactive UMETA(DisplayName = "Inactive Cell")
};

UCLASS()
class IVORYMEMORIES_API ACell : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ACell();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	//virtual void Tick(float DeltaTime) override;

	UPROPERTY(VisibleAnywhere, BlueprintReadWrite, Category = "Components")
	class UStaticMeshComponent* CellMesh; // Mesh for a cells

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cell")
	ECellState State; // Current cell state

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cell")
	FVector2D Coordinates; // X Y coords on the Board

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cell")
	AFigure* OccupiedFigure; // Reference to the piece occupying the cell (None if Epty)

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cell")
	FColor CellColor; // Cell color

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Cell")
	bool bIsHighlighted; // Highlighting cell possible state

	//------------------------------------------------------------

	UFUNCTION(BlueprintCallable, Category = "Cell")
	bool IsEmpty() const; // Check if cell is empty

	UFUNCTION(BlueprintCallable, Category = "Cell")
	void SetFigure(AFigure* newFigure); // Set a newer figure on empty cell 

	UFUNCTION(BlueprintCallable, Category = "Cell")
	void RemoveFigure(); // Remove figure from an occupied cell

	UFUNCTION(BlueprintCallable, Category = "Cell")
	void SetState(ECellState newState); // Change cell state

};
