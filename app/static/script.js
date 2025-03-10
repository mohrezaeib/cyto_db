document.addEventListener("DOMContentLoaded", function () {
    const selectAll = document.getElementById("selectAllFields");
    const fieldCheckboxes = document.querySelectorAll(".field-checkbox");
    const dropdownMenu = document.querySelector(".dropdown-menu");

    // Prevent dropdown from closing when clicking inside
    dropdownMenu.addEventListener("click", function (event) {
        event.stopPropagation();
    });

    // Handle "Select All" checkbox
    selectAll.addEventListener("change", function () {
        fieldCheckboxes.forEach(checkbox => {
            checkbox.checked = selectAll.checked;
        });
    });

    // Ensure "All Fields" updates dynamically
    function updateAllFields() {
        selectAll.checked = [...fieldCheckboxes].every(checkbox => checkbox.checked);
    }

    fieldCheckboxes.forEach(checkbox => {
        checkbox.addEventListener("change", function () {
            updateAllFields();
        });
    });

    // Run on load to set initial state
    updateAllFields();
});
