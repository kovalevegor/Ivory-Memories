// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Board.generated.h"

class ACell;

UCLASS()
class IVORYMEMORIES_API ABoard : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ABoard();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;
	void PostEditChangeProperty(FPropertyChangedEvent& FPropertyChangedEvent) override; // to be updated insede UE Editor

public:	
	// Called every frame
	//virtual void Tick(float DeltaTime) override;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Board")
	int32 Width;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Board")
	int32 Height;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Board")
	float CellSize; // cellsize for positioning

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Board")
	float CellSpacing; // space between cells

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Board")
	TSubclassOf<ACell> CellClass; // a class to spawn a cell

	UPROPERTY(VisibleAnywhere, BlueprintReadOnly, Category = "Board")
	TArray<ACell*> Cells;

	//------------------------------------------------------------

	UFUNCTION(BlueprintCallable, Category = "Board")
	void InitializeBoard(); // create the board

	UFUNCTION(BlueprintCallable, Category = "Board")
	ACell* GetCellAt(int32 X, int32 Y) const; // Get sell at coords XY

	UFUNCTION(BlueprintCallable, Category = "Board")
	void ClearBoard(); // Clear the board (destroy cells)


};
