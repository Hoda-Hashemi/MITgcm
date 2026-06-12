document.addEventListener("click", (event) => {
  const tab = event.target.closest("[data-tab-target]");
  if (tab) {
    const tabList = tab.closest(".tab-list");
    const tabs = tabList ? tabList.querySelectorAll("[data-tab-target]") : [];
    const tabsRoot = tab.closest(".field-tabs");
    const panels = tabsRoot ? tabsRoot.querySelectorAll("[data-tab-panel]") : [];
    const target = tab.getAttribute("data-tab-target");

    tabs.forEach((button) => {
      const active = button === tab;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-selected", String(active));
    });

    panels.forEach((panel) => {
      const active = panel.getAttribute("data-tab-panel") === target;
      panel.classList.toggle("is-active", active);
      panel.hidden = !active;
    });
    return;
  }

  const openButton = event.target.closest("[data-modal-src]");
  const modal = document.querySelector("[data-modal]");
  if (openButton && modal) {
    const image = modal.querySelector("[data-modal-img]");
    const caption = modal.querySelector("[data-modal-caption]");
    image.src = openButton.getAttribute("data-modal-src");
    image.alt = openButton.getAttribute("data-modal-caption") || "";
    caption.textContent = openButton.getAttribute("data-modal-caption") || "";
    modal.hidden = false;
    return;
  }

  if (event.target.closest("[data-modal-close]") || event.target === modal) {
    modal.hidden = true;
  }
});

document.addEventListener("keydown", (event) => {
  if (event.key !== "Escape") {
    return;
  }
  const modal = document.querySelector("[data-modal]");
  if (modal) {
    modal.hidden = true;
  }
});
