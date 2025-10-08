// Fill out your copyright notice in the Description page of Project Settings.


#include "CombatSystem/Figure.h"
#include "CombatSystem/Cell.h"

// Sets default values
AFigure::AFigure()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	//PrimaryActorTick.bCanEverTick = true;

	CurrentCell = nullptr;
	Rank = 0;
	MaxMoveDistance = 1;

	FigureMesh = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("FigureMesh"));
	RootComponent = FigureMesh;

	ApplyMeshSettings(); // Default mesh applying (if set in editor)
}

// Called when the game starts or when spawned
void AFigure::BeginPlay()
{
	Super::BeginPlay();
	ApplyMeshSettings();
	// Additional init if needed
	if (CurrentCell)
	{
		SetActorLocation(CurrentCell->GetActorLocation() + FVector(0.0f, 0.0f, 50.0f)); // Elevate the figure above the cell
	}
}

bool AFigure::CanMoveTo(ACell* TargetCell) const
{
	if (!TargetCell || TargetCell->State == ECellState::Inactive || !TargetCell->IsEmpty()) {
		return false; // occupied or inactive cell is not allowed
	}

	if (CurrentCell) {
		// Check the distance (Manhattan for a chess-like board)
		FVector2D Delta = TargetCell->Coordinates - CurrentCell->Coordinates;
		int32 Distance = FMath::Abs(Delta.X) + FMath::Abs(Delta.Y);
		return Distance <= MaxMoveDistance && Distance > 0;
	}
	return true; // if there is no current cell (default placement)
}

void AFigure::MoveToCell(ACell* newCell)
{
	if (CanMoveTo(newCell)) {
		if (CurrentCell) {
			CurrentCell->RemoveFigure();
		}
		SetCell(newCell);
		newCell->SetFigure(this);
		SetActorLocation(newCell->GetActorLocation() + FVector(0.0f, 0.0f, 50.0f));
	}
}

bool AFigure::CanAttack(AFigure* TargetFigure) const
{
	// By default any figure can attack any other figure if it exists and not itself
	return TargetFigure && TargetFigure != this;
}

void AFigure::Attack(AFigure* TargetFigure)
{
	if (CanAttack(TargetFigure)) {
		TargetFigure->Destroy();
		UE_LOG(LogTemp, Warning, TEXT("Figure %s attacked %s"), *FigureName, *TargetFigure->FigureName);
	}
}

//void AFigure::SetCell(ACell* newCell)
//{
//	CurrentCell = newCell;
//	FigureName = TEXT("Figure"); // base name must be edited by children
//	Rank = 0;
//	MaxMoveDistance = 1;
//	FigureMesh = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("FigureMesh"));
//	RootComponent = FigureMesh;
//	// HERE can be set default mesh, but have to be overrided by children
//}

void AFigure::SetCell(ACell* newCell)
{
	CurrentCell = newCell;
}

void AFigure::GenerateUniqueName(const FString& BaseName, int32 Index)
{
	FigureName = FString::Printf(TEXT("%s_%d"), *BaseName, Index);
}

void AFigure::ApplyMeshSettings()
{
	if (FigureMesh) {
		if (FigureMeshAsset) {
			FigureMesh->SetStaticMesh(FigureMeshAsset);
			UE_LOG(LogTemp, Warning, TEXT("Mesh applied: %s"), *FigureMeshAsset->GetName());
		}
		else {
			FigureMesh->SetStaticMesh(nullptr);
			UE_LOG(LogTemp, Warning, TEXT("Mesh cleared"));
		}

		if (FigureMaterial) { FigureMesh->SetMaterial(0, FigureMaterial); }
	}
}

// Automatic editing in Editor
void AFigure::PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent)
{
	Super::PostEditChangeProperty(PropertyChangedEvent);

	FName PropertyName = (PropertyChangedEvent.Property != nullptr) ?
		PropertyChangedEvent.Property->GetFName() : NAME_None;

	if (PropertyName == GET_MEMBER_NAME_CHECKED(AFigure, FigureMeshAsset) ||
		PropertyName == GET_MEMBER_NAME_CHECKED(AFigure, FigureMaterial)) {
		ApplyMeshSettings();
	}
}

//
//// Called every frame
//void AFigure::Tick(float DeltaTime)
//{
//	Super::Tick(DeltaTime);
//
//}

