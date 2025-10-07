// Fill out your copyright notice in the Description page of Project Settings.


#include "CombatSystem/Board.h"
#include "CombatSystem/Cell.h"
#include "Kismet/GameplayStatics.h" // ability to spawn

// Sets default values
ABoard::ABoard()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	//PrimaryActorTick.bCanEverTick = true;

	Width = 3;
	Height = 3;
	CellSize = 100.0f;
    CellSpacing = 10.0f;
	Cells.Empty();
}

// Called when the game starts or when spawned
void ABoard::BeginPlay()
{
	Super::BeginPlay();
	InitializeBoard(); // init when start the level
	
}

void ABoard::PostEditChangeProperty(FPropertyChangedEvent& PropertyChangedEvent)
{
	Super::PostEditChangeProperty(PropertyChangedEvent);
	// update the board in Editor when editing properties
	ClearBoard();
	InitializeBoard();
} 

void ABoard::InitializeBoard()
{
    if (!CellClass) return;
    UWorld* World = GetWorld();
    if (!World) return;

    // Инициализация массива
    Cells.SetNum(Width * Height);

    float BoardWidth = Width * (CellSize + CellSpacing) - CellSpacing;
    float BoardHeight = Height * (CellSize + CellSpacing) - CellSpacing;

    for (int32 Y = 0; Y < Height; Y++)
    {
        for (int32 X = 0; X < Width; X++)
        {
            float OffsetX = X * (CellSize + CellSpacing) - BoardWidth / 2.0f;
            float OffsetY = Y * (CellSize + CellSpacing) - BoardHeight / 2.0f;
            FVector Location = GetActorLocation() + FVector(OffsetX, OffsetY, 0.0f);
            FRotator Rotation = FRotator::ZeroRotator;
            FActorSpawnParameters SpawnParams;
            SpawnParams.Owner = this;

            ACell* NewCell = World->SpawnActor<ACell>(CellClass, Location, Rotation, SpawnParams);
            if (NewCell)
            {
                NewCell->Coordinates = FVector2D(X, Y);
                // Шахматная раскраска
                NewCell->CellColor = ((X + Y) % 2 == 0) ? FColor::White : FColor::Black;

                // Установка масштаба клетки на основе CellSize (предполагаем базовый меш 100x100)
                if (NewCell->CellMesh)
                {
                    float Scale = CellSize / 100.0f; // Если меш Plane — базовый размер 100
                    NewCell->CellMesh->SetRelativeScale3D(FVector(Scale, Scale, 1.0f));
                }

                // Прикрепление к доске (опционально)
                NewCell->AttachToActor(this, FAttachmentTransformRules::KeepWorldTransform);

                // Сохранение в массив
                int32 Index = Y * Width + X;
                Cells[Index] = NewCell;
            }
        }
    }
}

ACell* ABoard::GetCellAt(int32 X, int32 Y) const
{
    if (X >= 0 && X < Width && Y >= 0 && Y < Height)
    {
        int32 Index = Y * Width + X;
        return Cells[Index];
    }
    return nullptr;
}

void ABoard::ClearBoard()
{
    for (ACell* Cell : Cells)
    {
        if (Cell)
        {
            Cell->Destroy();
        }
    }
    Cells.Empty();
}

// Called every frame
//void ABoard::Tick(float DeltaTime)
//{
//	Super::Tick(DeltaTime);
//
//}

